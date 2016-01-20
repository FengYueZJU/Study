#include "RBEC.h"
#include "math.h"

double Potential::value(const double * p) const
{
    double result =  0.5 * (gamma_x * gamma_x * p[0] * p[0] + gamma_y * gamma_y * p[1] * p[1]);
    return result;
}

std::vector<double> Potential::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};


/// without normalization.

double Initial_Re::value(const double * p) const
{
    double result = ((1 - omega) + omega *  p[0]) / sqrt(PI) * exp(- (p[0] * p[0] + p[1] * p[1])/2);
    return result;
}

std::vector<double> Initial_Re::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};

double Initial_Im::value(const double * p) const
{
    double result = omega * p[1] / sqrt(PI) * exp(-(p[0] * p[0] + p[1] * p[1])/2);
    return result;
}

std::vector<double> Initial_Im::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};



RBEC::RBEC(const std::string& file) :
    mesh_file(file), beta(100.0), t(0.0), dt(1.0e-3), gamma_x(1.0), gamma_y(1.0), omega(0.75)
{};

RBEC::~RBEC()
{};

void RBEC::initialize()
{
    mesh.readData(mesh_file);

    template_geometry.readData("triangle.tmp_geo");
    coord_transform.readData("triangle.crd_trs");
    template_dof.reinit(template_geometry);
    template_dof.readData("triangle.2.tmp_dof");
    basis_function.reinit(template_dof);
    basis_function.readData("triangle.2.bas_fun");
    template_element.resize(1);
    template_element[0].reinit(template_geometry,
			       template_dof,
			       coord_transform,
			       basis_function);

    fem_space.reinit(mesh, template_element);
    int n_element = mesh.n_geometry(DIM);
    fem_space.element().resize(n_element);
    for (int i = 0;i < n_element;i ++)
	fem_space.element(i).reinit(fem_space, i, 0);
    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    phi_re.reinit(fem_space);
    phi_im.reinit(fem_space);
};

void RBEC::run()
{
    initialize();

    buildMatrixStruct();

    initialValue();

    FEMFunction <double, DIM> phi2(fem_space);

    mu = 0;

    do {
       
    	stepForward();

    	for (int i = 0; i < fem_space.n_dof(); ++i)
    	    phi2(i) = phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i);
    	phi2.writeOpenDXData("phi2.dx");

    	std::cout << "t  = " << t << std::endl;
    } while (t < 3);
};

void RBEC::initialValue()
{

    Initial_Re  phi_re_0(omega);
    Initial_Im  phi_im_0(omega);

    FEMFunction <double, DIM> phi_star_0(fem_space);

    Operator::L2Project(phi_re_0, phi_re, Operator::LOCAL_LEAST_SQUARE, 5);
    Operator::L2Project(phi_im_0, phi_im, Operator::LOCAL_LEAST_SQUARE, 5);

    int n_dof = fem_space.n_dof();

    for (int i = 0; i < n_dof; ++i)
	phi_star_0(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));
    
    double L2Phi_0 = Functional::L2Norm(phi_star_0, 10);

    std::cout << "L2Norm = " << L2Phi_0 << std::endl;

    for (int i = 0; i < n_dof; ++i)
    {
	phi_re(i) /= L2Phi_0;
	phi_im(i) /= L2Phi_0; 
    }

    for (int i = 0; i < n_dof; ++i)
	phi_star_0(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));
    
    L2Phi_0 = Functional::L2Norm(phi_star_0, 10);

    std::cout << "L2Norm = " << L2Phi_0 << std::endl;

    FEMFunction<double, DIM> vh(fem_space);
    Potential V(gamma_x, gamma_y);
    Operator::L2Project(V, vh, Operator::LOCAL_LEAST_SQUARE, 5);
    vh.writeOpenDXData("V.dx");
    phi_star_0.writeOpenDXData("phi0.dx");
};


void RBEC::stepForward()
{
    int i, j, k, l;
    int n_dof = fem_space.n_dof();
    int n_total_dof = 2 * n_dof;

    std::cout << "mu_old = " << mu << std::endl;

    mat_RBEC.reinit(sp_RBEC);
    mat_rere.reinit(sp_rere);
    mat_reim.reinit(sp_reim);
    mat_imre.reinit(sp_imre);
    mat_imim.reinit(sp_imim);

    Vector<double> phi(n_total_dof);
    FEMFunction <double, DIM> phi_star(fem_space);
    Vector<double> rhs(n_total_dof);
    Potential V(gamma_x, gamma_y);

/// 准备一个遍历全部单元的迭代器.
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();

/// 循环遍历全部单元, 只是为了统计每一行的非零元个数.
    for (; the_element != end_element; ++the_element)
    {
/// 当前单元信息.
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(10);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
/// 单元信息.
	std::vector<std::vector<std::vector<double> > > basis_gradient = the_element->basis_function_gradient(q_point);
	std::vector<std::vector<double> >  basis_value = the_element->basis_function_value(q_point);
	std::vector<double> phi_re_value = phi_re.value(q_point, *the_element);
	std::vector<double> phi_im_value = phi_im.value(q_point, *the_element);
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
/// 实际拼装.
	for (l = 0; l < n_quadrature_point; ++l)
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    for (j = 0; j < n_element_dof; ++j)
	    {
		for (k = 0; k < n_element_dof; ++k)
		{

		    double cont = Jxw * ((1 / dt) * basis_value[j][l] * basis_value[k][l]
					 + 0.5 * innerProduct(basis_gradient[j][l], basis_gradient[k][l])
					 + V.value(q_point[l]) * basis_value[j][l] * basis_value[k][l]
					 + beta * (phi_re_value[l] * phi_re_value[l]  + phi_im_value[l] * phi_im_value[l]) * basis_value[j][l] * basis_value[k][l] - mu * basis_value[j][l] * basis_value[k][l]);

		    mat_RBEC.add(element_dof[j], element_dof[k], cont);
		    mat_RBEC.add(element_dof[j] + n_dof, element_dof[k] + n_dof, cont);


		    cont = Jxw * (-omega * (q_point[l][0] * basis_gradient[k][l][1] 
					    - q_point[l][1] * basis_gradient[k][l][0]) * basis_value[j][l]);

		    mat_RBEC.add(element_dof[j], element_dof[k] + n_dof, cont);
		    mat_RBEC.add(element_dof[j] + n_dof, element_dof[k], -cont);
		}
		rhs(element_dof[j]) += Jxw * phi_re_value[l] * basis_value[j][l] / dt;
		rhs(element_dof[j] + n_dof) += Jxw * phi_im_value[l] * basis_value[j][l] / dt;
	    }
	}
    }

    
     FEMFunction<double, DIM> _phi_re(phi_re);
     FEMFunction<double, DIM> _phi_im(phi_im);

     boundaryValue(phi, rhs, mat_RBEC);


     dealii::SolverControl solver_control(4000, 1e-15);
     SolverGMRES<Vector<double> >::AdditionalData para(500, false, true);
     SolverGMRES<Vector<double> > gmres(solver_control, para);
     gmres.solve(mat_RBEC, phi, rhs, PreconditionIdentity());
        
     for (int i = 0; i < n_dof; ++i)
     {
        phi_re(i) = phi(i);
        phi_im(i) = phi(n_dof + i);
      }
   	

     for (int i = 0; i < n_dof; ++i)
	 phi_star(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));

     double L2Phi = Functional::L2Norm(phi_star, 10);

     std::cout << "L2 norm = " << L2Phi << std::endl;

     mu = -1/dt * log(L2Phi);
     std::cout << "mu_new = " << mu << std::endl;
    
     double e = energy(phi_re, phi_im, 10);
     std::cout << "Energy = " << e << std::endl;


     t += dt;
};

double RBEC::energy(FEMFunction<double, DIM>& phi_re, FEMFunction<double, DIM>& phi_im, int algebric_accuracy)
{
    double e = 0;
    Potential V(gamma_x, gamma_y);

    FEMSpace<double, DIM>& fem_space = phi_re.femSpace();
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) 
    {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<double> phi_re_value = phi_re.value(q_point, *the_element);
	std::vector<double> phi_im_value = phi_im.value(q_point, *the_element);
	std::vector<std::vector<double> > phi_re_gradient = phi_re.gradient(q_point, *the_element);
	std::vector<std::vector<double> > phi_im_gradient = phi_im.gradient(q_point, *the_element); 
	for (int l = 0; l < n_quadrature_point; l++) 
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    double phi2 = phi_re_value[l] * phi_re_value[l] + phi_im_value[l] * phi_im_value[l];

	    e += Jxw * (0.5 * (phi_re_gradient[l][0] * phi_re_gradient[l][0] +phi_re_gradient[l][1] * phi_re_gradient[l][1] + phi_im_gradient[l][0] * phi_im_gradient[l][0] +  phi_im_gradient[l][1] * phi_im_gradient[l][1]) 
			+ V.value(q_point[l]) * phi2
			+ 0.5 * beta * phi2 * phi2 - omega * (q_point[l][0] * phi_re_value[l] * phi_im_gradient[l][1]                
							      - q_point[l][1] * phi_re_value[l] * phi_im_gradient[l][0] 
							      - q_point[l][0] * phi_im_value[l] * phi_re_gradient[l][1] 
							      + q_point[l][1] * phi_im_value[l] * phi_im_gradient[l][0]));
	}
    }
    return e;
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////

