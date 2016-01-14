#include "BEC.h"

double Potential::value(const double * p) const
{
    double result =  0.5 * (gamma_x * gamma_x * p[0] * p[0] + gamma_y * gamma_y * p[1] * p[1]) + 4.0 * exp(-((p[0] - 1.0) * (p[0] - 1.0) + p[1] * p[1]));
    return result;
}

std::vector<double> Potential::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};

double u(const double * p)
{
    return 0;
};


/// without normalization.

double Initial::value(const double * p) const
{
    double result = sqrt(sqrt((gamma_x * gamma_y))) / sqrt(PI) * exp(- (p[0] * p[0] + p[1] * p[1])/2);
    return result;
}

std::vector<double> Initial_Re::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};


void RBEC::Matrix_A::getElementMatrix(
    const Element<double, DIM>& element0,
    const Element<double, DIM>& element1,
    const ActiveElementPairIterator<DIM>::State state)
{
    int n_element_dof0 = elementDof0().size();
    int n_element_dof1 = elementDof1().size();
    double volume = element0.templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebricAccuracy());
    std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double> > basis_value = element0.basis_function_value(q_point);
    std::vector<std::vector<std::vector<double> > > basis_gradient = element0.basis_function_gradient(q_point);
    std::vector<double>  phi_re_value = phi_re->value(q_point, element0);   
    std::vector<double>  phi_im_value = phi_im->value(q_point, element0); 
    Potential V(gamma_x, gamma_y);
    for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (int j = 0;j < n_element_dof0;j ++) {
	    for (int k = 0;k < n_element_dof1;k ++) {
		elementMatrix(j, k) += Jxw * ((1 / dt) * basis_value[j][l] * basis_value[k][l]
					      + 0.5 * innerProduct(basis_gradient[j][l], basis_gradient[k][l])
					      + V.value(q_point[l]) * basis_value[j][l] * basis_value[k][l]
					      + beta * (phi_re_value[l] * phi_re_value[l]  + phi_im_value[l] * phi_im_value[l]) * basis_value[j][l] * basis_value[k][l]);
	    }
	}
    }
};


BEC::BEC(const std::string& file) :
    mesh_file(file), beta(200.0), t(0.0), dt(1.0e-3), gamma_x(1.0), gamma_y(1.0)
{};

BEC::~BEC()
{};

void BEC::initialize()
{
    mesh.readData(mesh_file);

    template_geometry.readData("triangle.tmp_geo");
    coord_transform.readData("triangle.crd_trs");
    template_dof.reinit(template_geometry);
    template_dof.readData("triangle.1.tmp_dof");
    basis_function.reinit(template_dof);
    basis_function.readData("triangle.1.bas_fun");
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

    phi.reinit(fem_space);

};

void BEC::run()
{
    initialize();
    initialValue();
    phi.writeOpenDXData("phi0.dx");

    FEMFunction <double, DIM> phi2(fem_space);

    do {
	stepForward();
	phi.writeOpenDXData("phi.dx");

	for (int i = 0; i < fem_space.n_dof(); ++i)
	    phi2(i) = phi(i) * phi(i);
	phi2.writeOpenDXData("phi2.dx");

	std::cout << "t  = " << t << std::endl;
    } while (t < 3);
};

void RBEC::initialValue()
{
    Initial  phi_0(gamma_x,gamma_y);

//    FEMFunction <double, DIM> phi_star_0(fem_space);

    Operator::L2Project(phi_0, phi, Operator::LOCAL_LEAST_SQUARE, 3);

    int n_dof = fem_space.n_dof();
    ///// the same as .* in matlab
    for (int i = 0; i < n_dof; ++i)
	phi_star_0(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));
    
    double L2Phi_0 = Functional::L2Norm(phi_star_0, 6);
    std::cout << "L2Norm = " << L2Phi_0 << std::endl;

    for (int i = 0; i < n_dof; ++i)
    {
	phi_re(i) /= L2Phi_0;
	phi_im(i) /= L2Phi_0; 
    }

    FEMFunction<double, DIM> vh(fem_space);
    Potential V(gamma_x, gamma_y);
    Operator::L2Project(V, vh, Operator::LOCAL_LEAST_SQUARE, 3);
    vh.writeOpenDXData("V.dx");
};


void RBEC::stepForward()
{
    int i, j, k, l;
    int n_dof_phi = fem_space.n_dof();

    Matrix_A matA(fem_space, dt, beta, gamma_x, gamma_y, omega, phi_re , phi_im);
    Matrix_B matB(fem_space, dt, beta, gamma_x, gamma_y, omega, phi_re , phi_im);
    // Matrix_C matC(fem_space, dt, beta, gamma_x, gamma_y, omega, phi_re , phi_im);
    Matrix_D matD(fem_space, dt, beta, gamma_x, gamma_y, omega, phi_re , phi_im);

    matA.algebricAccuracy() = 6;
    matA.build();
    matB.algebricAccuracy() = 6;
    matB.build();
    matD.algebricAccuracy() = 6;
    matD.build();

    Vector<double> rhs_re(n_dof_phi);
    Vector<double> rhs_im(n_dof_phi);

    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;++ the_element) {
    	double volume = the_element->templateElement().volume();
    	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(DIM);
    	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    	int n_quadrature_point = quad_info.n_quadraturePoint();
    	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    	std::vector<double> phi_re_value = phi_re.value(q_point, *the_element);
    	std::vector<double> phi_im_value = phi_im.value(q_point, *the_element);
    	int n_element_dof = the_element->n_dof();
    	const std::vector<int>& element_dof = the_element->dof();
    	for (l = 0;l < n_quadrature_point;l ++) 
    	{
    	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
    	    for (j = 0;j < n_element_dof;j ++) 
    	    {
    		rhs_re(element_dof[j]) += Jxw * phi_re_value[l] * basis_value[j][l] / dt;
    		rhs_im(element_dof[j]) += Jxw * phi_im_value[l] * basis_value[j][l] / dt;
    	    }
    	}
    }
    double err = 1.0;
    while (err > 1e-10)
    {
	FEMFunction<double, DIM> _phi_re(phi_re);
	FEMFunction<double, DIM> _phi_im(phi_im);
	Vector<double> tmp1(phi_im);
	matB.vmult(tmp1, tmp1);
    
	for (int i = 0; i < n_dof_phi; ++i)
	    rhs_re(i) -= tmp1(i);

	BoundaryFunction<double, DIM> boundary_A(BoundaryConditionInfo::DIRICHLET, 1, &u_re);
	BoundaryConditionAdmin<double, DIM> boundary_admin_A(fem_space);
	boundary_admin_A.add(boundary_A);
	boundary_admin_A.apply(matA, phi_re, rhs_re);

	AMGSolver solverA(matA);
	solverA.solve(phi_re, rhs_re);

	Vector<double> tmp2(phi_re);
	matB.vmult(tmp2, tmp2);
    
	for (int i = 0; i < n_dof_phi; ++i)
	    rhs_im(i) += tmp2(i);
    
	BoundaryFunction<double, DIM> boundary_B(BoundaryConditionInfo::DIRICHLET, 1, &u_im);
	BoundaryConditionAdmin<double, DIM> boundary_admin_B(fem_space);
	boundary_admin_B.add(boundary_B);
	boundary_admin_B.apply(matD, phi_im, rhs_im);

//   AMGSolver solverA(matB);
	solverA.solve(phi_im, rhs_im);

//    FEMFunction<double, DIM> delta_phi_re(fem_space);
//    FEMFunction<double, DIM> delta_phi_im(fem_space);
//	FEMFunction<double, DIM> phi_star(fem_space);
	err = 0;

	for (int i = 0; i < n_dof_phi; ++i)
	{
	    err += (_phi_re(i) - phi_re(i)) * (_phi_re(i) - phi_re(i));  
	    err += (_phi_im(i) - phi_im(i)) * (_phi_im(i) - phi_im(i));
	}
	err = sqrt(err);
  
	std::cout << "err = " << err << std::endl;
    // for (int i = 0; i < n_dof_phi; ++i)
    // 	phi_star(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));

    }
    // for (int i = 0; i < n_dof_phi; ++i)
    // 	phi_star(i) = sqrt(phi_re(i) * phi_re(i) + phi_im(i) * phi_im(i));

    double L2Phi = Functional::L2Norm(phi_re, 6);

    std::cout << "L2 norm = " << L2Phi << std::endl;
    for (int i = 0; i < n_dof_phi; ++i)
    {
    	phi_re(i) /= L2Phi;
    	phi_im(i) /= L2Phi;
//    	delta_phi_re(i) -= phi_re(i);
//    	delta_phi_im(i) -= phi_im(i);
//    	delta_ phi_star(i) = sqrt(delta_phi_re(i) * delta_phi_re(i) + delta_phi_im(i) * delta_phi_im(i)); 
    }

//    double res_phi = Functional::L2Norm(delta_phi_star, 6);

//    double e = energy(phi_re,phi_im, 6);
//    std::cout << "Energy = " << e << std::endl;
//    std::cout << "Res = " << res_phi << std::endl;
    
    t += dt;
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
//////////Energy of Lz
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

