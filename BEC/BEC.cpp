#include "BEC.h"

double Potential::value(const double * p) const
{
    double result =  0.5 * (gamma_x * gamma_x * p[0] * p[0] + gamma_y * gamma_y * p[1] * p[1]) + omega0 * exp(-delta *((p[0] - r0) * (p[0] - r0) + p[1] * p[1]));
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

double Initial::value(const double * p) const
{
    double result = pow((gamma_x * gamma_y),0.25) / sqrt(PI) * exp(- (gamma_x * p[0] * p[0] + gamma_y *  p[1] * p[1])/2);
    return result;
}

std::vector<double> Initial::gradient(const double * p) const
{
    std::vector<double> v(DIM);
    return v;
};


void BEC::Matrix::getElementMatrix(
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
    std::vector<double>  phi_value = phi->value(q_point, element0);   
    Potential V(gamma_x, gamma_y,omega0,delta,r0);
    for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (int j = 0;j < n_element_dof0;j ++) {
	    for (int k = 0;k < n_element_dof1;k ++) {
		elementMatrix(j, k) += Jxw * ((1 / dt) * basis_value[j][l] * basis_value[k][l]
					      + 0.5 * innerProduct(basis_gradient[j][l], basis_gradient[k][l])
					      + V.value(q_point[l]) * basis_value[j][l] * basis_value[k][l]
					      + beta * (phi_value[l] * phi_value[l]) * basis_value[j][l] * basis_value[k][l]);
	    }
	}
    }
};


BEC::BEC(const std::string& file) :
    mesh_file(file), beta(200.0), t(0.0), dt(1.0e-3), gamma_x(1.0), gamma_y(1.0), omega0(4.0), delta(1.0), r0(1.0)
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

void BEC::initialValue()
{
    Initial  phi_0(gamma_x,gamma_y);

    Operator::L2Project(phi_0, phi, Operator::LOCAL_LEAST_SQUARE, 3);

    int n_dof = fem_space.n_dof();
    
    double L2Phi_0 = Functional::L2Norm(phi, 6);
    std::cout << "L2Norm = " << L2Phi_0 << std::endl;

    for (int i = 0; i < n_dof; ++i)
    {
	phi(i) /= L2Phi_0; 
    }

    FEMFunction<double, DIM> vh(fem_space);
    Potential V(gamma_x, gamma_y,omega0,delta,r0);
    Operator::L2Project(V, vh, Operator::LOCAL_LEAST_SQUARE, 3);
    vh.writeOpenDXData("V.dx");
};


void BEC::stepForward()
{
    int i, j, k, l;
    int n_dof_phi = fem_space.n_dof();

    Matrix mat(fem_space, dt, beta, gamma_x, gamma_y, omega0, delta, r0, phi);

    mat.algebricAccuracy() = 6;
    mat.build();

    Vector<double> rhs(n_dof_phi);

    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;++ the_element) {
    	double volume = the_element->templateElement().volume();
    	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(DIM);
    	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    	int n_quadrature_point = quad_info.n_quadraturePoint();
    	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    	std::vector<double> phi_value = phi.value(q_point, *the_element);
    	int n_element_dof = the_element->n_dof();
    	const std::vector<int>& element_dof = the_element->dof();
    	for (l = 0;l < n_quadrature_point;l ++) 
    	{
    	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
    	    for (j = 0;j < n_element_dof;j ++) 
    	    {
    		rhs(element_dof[j]) += Jxw * phi_value[l] * basis_value[j][l] / dt;
    	    }
    	}
    }
    double err = 1.0;
    while (err > 1e-10)
    {
	FEMFunction<double, DIM> _phi(phi);

	BoundaryFunction<double, DIM> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
	BoundaryConditionAdmin<double, DIM> boundary_admin(fem_space);
	boundary_admin.add(boundary);
	boundary_admin.apply(mat, phi, rhs);

	AMGSolver solver(mat);
	solver.solve(phi, rhs);

	err = 0;

	for (int i = 0; i < n_dof_phi; ++i)
	{
	    err += (_phi(i) - phi(i)) * (_phi(i) - phi(i));  
	}
	err = sqrt(err);
  
	std::cout << "err = " << err << std::endl;
  
    }
        
    FEMFunction<double, DIM> delta_phi(fem_space);
 
    double L2Phi = Functional::L2Norm(phi, 6);

    std::cout << "L2 norm = " << L2Phi << std::endl;

    for (int i = 0; i < n_dof_phi; ++i)
    {
    	phi(i) /= L2Phi;
    	delta_phi(i) -= phi(i);
    }
    double res_phi = Functional::L2Norm(delta_phi, 6);

    double e = energy(phi,6);
    std::cout << "Energy = " << e << std::endl;
    std::cout << "Res = " << res_phi << std::endl;
    
    t += dt;
};

double BEC::energy(FEMFunction<double, DIM>& phi, int algebric_accuracy)
{
    double e = 0;
    Potential V(gamma_x, gamma_y, omega0, delta, r0);

    FEMSpace<double, DIM>& fem_space = phi.femSpace();
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) 
    {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<double> phi_value = phi.value(q_point, *the_element);
	std::vector<std::vector<double> > phi_gradient = phi.gradient(q_point, *the_element);

	for (int l = 0; l < n_quadrature_point; l++) 
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    double phi2 = phi_value[l] * phi_value[l];

	    e += Jxw * (0.5 * (phi_gradient[l][0] * phi_gradient[l][0] +phi_gradient[l][1] * phi_gradient[l][1]) 
			+ V.value(q_point[l]) * phi2
			+ 0.5 * beta * phi2 * phi2);

	}
    }
    return e;
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////

