#include "RBEC.h"

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
    double result = ((1.0 - omega) + omega *  p[0]) / sqrt(PI) * exp(- (p[0] * p[0] + p[1] * p[1])/2.0);
//    double result = 1.0 / sqrt(PI) *  exp(- (p[0] * p[0] + p[1] * p[1])/2.0); 
    return result;
}

std::vector<double> Initial_Re::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};

double Initial_Im::value(const double * p) const
{
    double result = omega * p[1] / sqrt(PI) * exp(-(p[0] * p[0] + p[1] * p[1])/2.0);
//    double result = 0;  
    return result;
}

std::vector<double> Initial_Im::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};



RBEC::RBEC(const std::string& file) :
    mesh_file(file), beta(100.0), t(0.0), dt(1.0e-3), gamma_x(1.0), gamma_y(1.0), omega(0.75), Linf_delta_phi(1.0)
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
    initialMesh();
};

void RBEC::run()
{
    initialize();

    buildMatrixStruct();

    do {

	moveMesh();
    	stepForward();
	outputSolution();

    	std::cout << "t  = " << t << std::endl;
    } while (Linf_delta_phi > 1e-7);
};


void RBEC::initialMesh()
{
    std::cout << "Initialize mesh ... " << std::endl;
    double scale, scale_step = 0.2;
    scale = scale_step;
    FEMFunction <double, DIM> phi2(fem_space);
    do {
	  initialValue();
	  phi2.scale(scale);
	  moveMesh();
	  outputSolution();
	  std::cout << "\r\tscale = " << scale << std::endl;
	  scale += scale_step;
	} while (scale <= 1.0);
};

void RBEC::initialValue()
{

    Initial_Re  phi_re_0(omega);
    Initial_Im  phi_im_0(omega);

    FEMFunction <double, DIM> phi_star_0(fem_space);
    FEMFunction <double, DIM> phi2(fem_space);
  
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
        phi2(i) = phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i);
	
    FEMFunction<double, DIM> vh(fem_space);
    Potential V(gamma_x, gamma_y);
    Operator::L2Project(V, vh, Operator::LOCAL_LEAST_SQUARE, 5);
    vh.writeOpenDXData("V.dx");
    phi2.writeOpenDXData("phi0.dx");
};


void RBEC::stepForward()
{
    int i, j, k, l;
    int n_dof = fem_space.n_dof();
    int n_total_dof = 2 * n_dof;

    mat_RBEC.reinit(sp_RBEC);
//    mat_rere.reinit(sp_rere);
//    mat_reim.reinit(sp_reim);
//    mat_imre.reinit(sp_imre);
//    mat_imim.reinit(sp_imim);

    Vector<double> phi(n_total_dof);
    FEMFunction <double, DIM> phi_star(fem_space);
    FEMFunction <double, DIM> _phi_re(phi_re);
    FEMFunction <double, DIM> _phi_im(phi_im);
    FEMFunction <double, DIM> delta_phi(fem_space);
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
					 + beta * (phi_re_value[l] * phi_re_value[l]  + phi_im_value[l] * phi_im_value[l]) * basis_value[j][l] * basis_value[k][l]);

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

    boundaryValue(phi, rhs, mat_RBEC);

//    const std::size_t * rowstart = sp_RBEC.get_rowstart_indices();
//    const unsigned int * colnum = sp_RBEC.get_column_numbers();

//    std::ofstream output("Matrix.m");
//    output.setf(std::ios::fixed);
//    output.precision(20);
//    output << "A =[" << std::endl;
//	for (int i = 0; i < n_total_dof; ++i)
//	{
//		for (int j = 0; j < n_total_dof; ++j)
//			output << mat_RBEC.el(i,j) << "\t";
//		output << std::endl;
//	}
//	output << "];" << std::endl;
//	output.close();
//	std::cout << "RBEC matrix outputed!" << std::endl;

//    std::ofstream output_rhs("rhs.m");
//    output_rhs.setf(std::ios::fixed);
//    output_rhs.precision(20);
//    output_rhs << "rhs =[" << std::endl;
//	for (int i = 0; i < n_total_dof; ++i)
//	{
//	    output_rhs << rhs(i) << std::endl;
//	}
//	output_rhs << "];" << std::endl;
//	output_rhs.close();
//	std::cout << "RBEC rhs outputed!" << std::endl;
//	getchar();
 

     PreconditionSSOR<> preconditioner;
     preconditioner.initialize(mat_RBEC, 1.2);

//     SparseILU<double> preconditioner;
//     SparseILU<> ilu;
//     ilu.initialize(mat_RBEC);

     dealii::SolverControl solver_control(200000, 1e-15);
//     SolverGMRES<Vector<double> >::AdditionalData para(500, false, true);
//     SolverGMRES<Vector<double> > gmres(solver_control, para);
     SolverGMRES<Vector<double> > gmres(solver_control);
     gmres.solve(mat_RBEC, phi, rhs, preconditioner);
//     gmres.solve(mat_RBEC, phi, rhs, PreconditionIdentity());
      
     for (int i = 0; i < n_dof; ++i)
     {
        phi_re(i) = phi(i);
        phi_im(i) = phi(n_dof + i);
      }
   	
     for (int i = 0; i < n_dof; ++i)
	 phi_star(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));

     double L2Phi = Functional::L2Norm(phi_star, 10);

     std::cout << "L2 norm = " << L2Phi << std::endl;
       
     for (int i = 0; i < n_dof; ++i)
     {
    	phi_re(i) /= L2Phi;
	phi_im(i) /= L2Phi;
	delta_phi(i) = sqrt((_phi_re(i) - phi_re(i)) * (_phi_re(i) - phi_re(i)) +  (_phi_im(i) - phi_im(i)) *  (_phi_im(i) - phi_im(i)));
     }

     Linf_delta_phi = delta_phi(0);

     for (int i = 1; i < n_dof; ++i)
     {
	 if (delta_phi(i) > Linf_delta_phi)
	     Linf_delta_phi = delta_phi(i);
     }
     std::cout << "Linf norm of delta_phi  = " << Linf_delta_phi << std::endl;

     double e = energy(phi_re, phi_im, 10);
     std::cout << "Energy = " << e << std::endl;

     t += dt;
};


void RBEC::getMonitor()
{
	int i, l;
        FEMFunction <double, DIM> phi2(fem_space);
	FEMSpace<double,DIM>::ElementIterator the_element = fem_space.beginElement();
	FEMSpace<double,DIM>::ElementIterator end_element = fem_space.endElement();
	for (i = 0;the_element != end_element;++ the_element) {
	        double volume = the_element->templateElement().volume();
           	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(10);
	        std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	        int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
		std::vector<std::vector<double> > phi2_gradient = phi2.gradient(q_point, *the_element); 

		float d = 0, area = 0;
		for (l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			area += Jxw;
			d += Jxw*innerProduct(phi2_gradient[l], phi2_gradient[l]);
		}
//	        std::cout << "d=" << d << std::endl;
//	        std::cout << "area=" << area << std::endl;
//		std::cout << "i = " << i << std::endl;
		monitor(i ++) = d/area;
	
	}
	std::cout << "max monitor=" << *std::max_element(monitor().begin(), monitor().end())
		<< "\tmin monitor=" << *std::min_element(monitor().begin(), monitor().end())
		<< std::endl;
	smoothMonitor(DIM);
	for (i = 0;i < n_geometry(DIM);i ++)
		monitor(i) = 1./sqrt(1. + monitor(i));
};



void RBEC::updateSolution()
{
	int i, j, l;
        int n_dof = fem_space.n_dof();
	int n_total_dof = 2 * n_dof;
	Vector<double> rhs(n_total_dof);
	Vector<double> phi(n_total_dof);
	FEMFunction<double,DIM> _phi_re(phi_re);
        FEMFunction<double,DIM> _phi_im(phi_im);
	const double& msl = moveStepLength();
//////////////////!!!!
	MassMatrix<DIM,double> matrix(fem_space);
	matrix.algebricAccuracy() = 2;
	matrix.build();
     
        PreconditionSSOR<> preconditioner;
        preconditioner.initialize(mat_RBEC, 1.2);

        dealii::SolverControl solver_control(200000, 1e-15);
        SolverGMRES<Vector<double> > gmres(solver_control);
        gmres.solve(mat_RBEC, phi, rhs, preconditioner);

	for (i = 1;i > 0;i --) {  

		FEMSpace<double,2>::ElementIterator the_element = fem_space.beginElement();
		FEMSpace<double,2>::ElementIterator end_element = fem_space.endElement();
		for (;the_element != end_element;++ the_element) {
			double volume = the_element->templateElement().volume();
			const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(2);
			std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
			int n_quadrature_point = quad_info.n_quadraturePoint();
			std::vector<Point<2> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
			std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
			std::vector<double> _phi_re_value = _phi_re.value(q_point, *the_element);
                        std::vector<double> _phi_im_value = _phi_im.value(q_point, *the_element);
			std::vector<std::vector<double> > phi_re_gradient = phi_re.gradient(q_point, *the_element);
	                std::vector<std::vector<double> > phi_im_gradient = phi_im.gradient(q_point, *the_element);
			std::vector<std::vector<double> > move_vector = moveDirection(q_point, the_element->index());
			int n_element_dof = the_element->n_dof();
			const std::vector<int>& element_dof = the_element->dof();
			for (l = 0;l < n_quadrature_point;l ++) {
				double Jxw = quad_info.weight(l)*jacobian[l]*volume;
				for (j = 0;j < n_element_dof;j ++) {
					rhs(element_dof[j]) += Jxw*basis_value[j][l]*(_phi_re_value[l]
							+ (1./i)*msl*innerProduct(move_vector[l], phi_re_gradient[l]));
	                                rhs(element_dof[j] + n_dof) += Jxw*basis_value[j][l]*(_phi_im_value[l]
							+ (1./i)*msl*innerProduct(move_vector[l], phi_im_gradient[l]));
				}
			}
		}

		gmres.solve(mat_RBEC, phi, rhs, preconditioner);
		for (int s = 0; s < n_dof; ++s)
		{
		    phi_re(s) = phi(s);
		    phi_im(s) = phi(s + n_dof);
                }
	};
};

void RBEC::outputSolution()
{

    FEMFunction <double, DIM> phi2(fem_space);

    outputPhysicalMesh("E");
    for (int i = 0; i < fem_space.n_dof(); ++i)
       phi2(i) = phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i);
    phi2.writeOpenDXData("phi2.dx");
};

double RBEC::energy(FEMFunction<double, DIM>& phi_re, FEMFunction<double, DIM>& phi_im, int algebric_accuracy)
{
    double e = 0;
    Potential V(gamma_x, gamma_y);
//////
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
							      + q_point[l][1] * phi_im_value[l] * phi_re_gradient[l][0]));

	}
    }
    return e;
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////

