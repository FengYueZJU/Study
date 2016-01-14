/**
 * @file   RBEC.h
 * @author Feng Yue <fengyue@zju.edu.cn>
 * @date   Wed Dec 23 09:58:53 2015
 * 
 * @brief  带旋转的BEC方程的类型声明.
 * 
 * 
 */

#include <AFEPack/MovingMesh2D.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/BilinearOperator.h>

#define DIM 2
#define PI (4.0 * atan(1.0))

/**
 * 势能函数类型.
 * 
 */
class Potential : public Function<double>
{
private:
    double gamma_x;		/**< 参数. */
    double gamma_y;		/**< 参数. */
public:
    /** 
     * 初始构造函数, 传递参数.
     * 
     * @param _gamma_x 参数. 
     * @param _gamma_y 参数.
     * 
     * @return 标准返回.
     */
    Potential(double _gamma_x, 
	      double _gamma_y) :
	gamma_x(_gamma_x), 
	gamma_y(_gamma_y)
	{};
    /** 
     * 析构函数.
     * 
     * 
     * @return 
     */
    ~Potential() {};
public:
    /** 
     * 求函数值.
     * 
     * 
     * @return 函数值. 
     */
    virtual double value(const double *) const;

    /** 
     * 求梯度值.
     * 
     * 
     * @return 梯度值.
     */
    virtual std::vector<double> gradient(const double *) const;
};

/**
 *  初值函数类型, 实部.
 * 
 */
class Initial_Re : public Function<double>
{
private:
    double omega;		/**< 角速度. */
public:
    /** 
     * 初始构造函数.
     * 
     * @param _omega 角速度.
     * 
     * @return 标准返回.
     */
    Initial_Re(double _omega) :
	omega(_omega)
	{};
    /** 
     * 析构函数.
     * 
     * 
     * @return 标准返回. 
     */
    ~Initial_Re() {};
public:
    /** 
     * 求函数值.
     * 
     * 
     * @return 函数值.
     */
    virtual double value(const double *) const;
    /** 
     * 求梯度值.
     * 
     * 
     * @return 梯度值.
     */
    virtual std::vector<double> gradient(const double *) const;
};

/**
 * 初值函数类型, 虚部.
 * 
 */
class Initial_Im : public Function<double>
{
private:
    double omega;		/**< 角速度. */
public:

    /** 
     * 初始构造函数.
     * 
     * @param _omega 角速度. 
     * 
     * @return 标准返回.
     */
    Initial_Im(double _omega) :
	omega(_omega)
	{};
    /** 
     * 析构函数.
     * 
     * 
     * @return 标准返回.
     */
    ~Initial_Im() {};
public:
    /** 
     * 求函数值.
     * 
     * 
     * @return 函数值. 
     */
    virtual double value(const double *) const;
    /** 
     * 求梯度值.
     * 
     * 
     * @return 梯度值.
     */
    virtual std::vector<double> gradient(const double *) const;
};

/**
 * 带旋转的BEC方程类.
 * 
 */
class RBEC
{
public:
    /**
     * 系数矩阵类. 解空间实部, 测试空间实部.
     * 
     */
    class Matrix_A : public StiffMatrix<DIM, double>
    {
    private:
	double dt;		/**< 时间步长. */
	double beta;		/**< 参数. */
	double omega;		/**< 角速度. */
	double gamma_x;		/**< 参数. */
	double gamma_y;		/**< 参数. */
	FEMFunction<double, DIM> *phi_re; /**< 数值解实部. */
	FEMFunction<double, DIM> *phi_im; /**< 数值解虚部. */
    public:
	/** 
	 * 构造函数.
	 * 
	 * @param sp 有限元空间.
	 * @param _dt 时间步长.
	 * @param _beta 参数.
	 * @param _omega 角速度.
	 * @param _gamma_x 参数.
	 * @param _gamma_y 参数.
	 * @param _phi_re 数值解, 实部.
	 * @param _phi_im 数值解, 虚部.
	 * 
	 * @return 标准返回.
	 */
	Matrix_A(FEMSpace<double, DIM>& sp, 
		 const double& _dt, 
		 const double& _beta,
		 const double& _omega,
		 const double& _gamma_x,
		 const double& _gamma_y,
		 FEMFunction<double, DIM>& _phi_re,
		 FEMFunction<double, DIM>& _phi_im) :
	    dt(_dt), 
	    beta(_beta),	
            omega(_omega),
	    gamma_x(_gamma_x),
	    gamma_y(_gamma_y),
	    StiffMatrix<DIM, double>(sp) 
	    {
		phi_re = &_phi_re;
		phi_im = &_phi_im;
	    };
	/** 
	 * 析构函数.
	 * 
	 * 
	 * @return 标准返回. 
	 */
	virtual ~Matrix_A() {};
    public:
	/** 
	 * 拼装矩阵元素.
	 * 
	 * @param e0 解空间单元.
	 * @param e1 测试空间单元.
	 * @param state 拼装状态.
	 */
	virtual void getElementMatrix(const Element<double, DIM>& e0,
				      const Element<double, DIM>& e1,
				      const ActiveElementPairIterator<DIM>::State state);
    };

    /**
     * 系数矩阵, 解空间实部, 测试空间虚部. 
     * 
     */
    class Matrix_B : public StiffMatrix<DIM, double>
    {
    private:
	double dt;		/**< 时间步长. */
	double beta;		/**< 参数. */
	double omega;		/**< 角速度. */
	double gamma_x;		/**< 参数. */
	double gamma_y;		/**< 参数. */
	FEMFunction<double, DIM> *phi_re; /**< 数值解实部. */
	FEMFunction<double, DIM> *phi_im; /**< 数值解虚部. */
    public:
	/** 
	 * 构造函数.
	 * 
	 * @param sp 有限元空间.
	 * @param _dt 时间步长.
	 * @param _beta 参数.
	 * @param _omega 角速度.
	 * @param _gamma_x 参数.
	 * @param _gamma_y 参数.
	 * @param _phi_re 数值解, 实部.
	 * @param _phi_im 数值解, 虚部.
	 * 
	 * @return 标准返回.
	 */
	Matrix_B(FEMSpace<double, DIM>& sp, 
		 const double& _dt, 
		 const double& _beta,
		 const double& _omega,
		 const double& _gamma_x,
		 const double& _gamma_y,
		 FEMFunction<double, DIM>& _phi_re,
		 FEMFunction<double, DIM>& _phi_im) :
	    dt(_dt), 
	    beta(_beta),	
            omega(_omega),
	    gamma_x(_gamma_x),
	    gamma_y(_gamma_y),
	    StiffMatrix<DIM, double>(sp) 
	    {
		phi_re = &_phi_re;
		phi_im = &_phi_im;
	    };
	/** 
	 * 析构函数.
	 * 
	 * 
	 * @return 标准返回. 
	 */
	virtual ~Matrix_B() {};
    public:
	/** 
	 * 拼装矩阵元素.
	 * 
	 * @param e0 解空间单元.
	 * @param e1 测试空间单元.
	 * @param state 拼装状态.
	 */
	virtual void getElementMatrix(const Element<double, DIM>& e0,
				      const Element<double, DIM>& e1,
				      const ActiveElementPairIterator<DIM>::State state);
    };

    /**
     * 系数矩阵类. 解空间实部, 测试空间实部.
     * 
     */
    class Matrix_C : public StiffMatrix<DIM, double>
    {
    private:
	double dt;		/**< 时间步长. */
	double beta;		/**< 参数. */
	double omega;		/**< 角速度. */
	double gamma_x;		/**< 参数. */
	double gamma_y;		/**< 参数. */
	FEMFunction<double, DIM> *phi_re; /**< 数值解实部. */
	FEMFunction<double, DIM> *phi_im; /**< 数值解虚部. */
    public:
	/** 
	 * 构造函数.
	 * 
	 * @param sp 有限元空间.
	 * @param _dt 时间步长.
	 * @param _beta 参数.
	 * @param _omega 角速度.
	 * @param _gamma_x 参数.
	 * @param _gamma_y 参数.
	 * @param _phi_re 数值解, 实部.
	 * @param _phi_im 数值解, 虚部.
	 * 
	 * @return 标准返回.
	 */
	Matrix_C(FEMSpace<double, DIM>& sp, 
		 const double& _dt, 
		 const double& _beta,
		 const double& _omega,
		 const double& _gamma_x,
		 const double& _gamma_y,
		 FEMFunction<double, DIM>& _phi_re,
		 FEMFunction<double, DIM>& _phi_im) :
	    dt(_dt), 
	    beta(_beta),	
            omega(_omega),
	    gamma_x(_gamma_x),
	    gamma_y(_gamma_y),
	    StiffMatrix<DIM, double>(sp) 
	    {
		phi_re = &_phi_re;
		phi_im = &_phi_im;
	    };
	/** 
	 * 析构函数.
	 * 
	 * 
	 * @return 标准返回. 
	 */
	virtual ~Matrix_C() {};
    public:
	/** 
	 * 拼装矩阵元素.
	 * 
	 * @param e0 解空间单元.
	 * @param e1 测试空间单元.
	 * @param state 拼装状态.
	 */
	virtual void getElementMatrix(const Element<double, DIM>& e0,
				      const Element<double, DIM>& e1,
				      const ActiveElementPairIterator<DIM>::State state);
    };
    /**
     * 系数矩阵类. 解空间实部, 测试空间实部.
     * 
     */
    class Matrix_D : public StiffMatrix<DIM, double>
    {
    private:
	double dt;		/**< 时间步长. */
	double beta;		/**< 参数. */
	double omega;		/**< 角速度. */
	double gamma_x;		/**< 参数. */
	double gamma_y;		/**< 参数. */
	FEMFunction<double, DIM> *phi_re; /**< 数值解实部. */
	FEMFunction<double, DIM> *phi_im; /**< 数值解虚部. */
    public:
	/** 
	 * 构造函数.
	 * 
	 * @param sp 有限元空间.
	 * @param _dt 时间步长.
	 * @param _beta 参数.
	 * @param _omega 角速度.
	 * @param _gamma_x 参数.
	 * @param _gamma_y 参数.
	 * @param _phi_re 数值解, 实部.
	 * @param _phi_im 数值解, 虚部.
	 * 
	 * @return 标准返回.
	 */
	Matrix_D(FEMSpace<double, DIM>& sp, 
		 const double& _dt, 
		 const double& _beta,
		 const double& _omega,
		 const double& _gamma_x,
		 const double& _gamma_y,
		 FEMFunction<double, DIM>& _phi_re,
		 FEMFunction<double, DIM>& _phi_im) :
	    dt(_dt), 
	    beta(_beta),	
            omega(_omega),
	    gamma_x(_gamma_x),
	    gamma_y(_gamma_y),
	    StiffMatrix<DIM, double>(sp) 
	    {
		phi_re = &_phi_re;
		phi_im = &_phi_im;
	    };
	/** 
	 * 析构函数.
	 * 
	 * 
	 * @return 标准返回. 
	 */
	virtual ~Matrix_D() {};
    public:
	/** 
	 * 拼装矩阵元素.
	 * 
	 * @param e0 解空间单元.
	 * @param e1 测试空间单元.
	 * @param state 拼装状态.
	 */
	virtual void getElementMatrix(const Element<double, DIM>& e0,
				      const Element<double, DIM>& e1,
				      const ActiveElementPairIterator<DIM>::State state);
    };

private:
    TemplateGeometry<DIM> template_geometry; /**< 参考几何体. */
    CoordTransform<DIM, DIM> coord_transform; /**< 座标变换. */
    TemplateDOF<DIM> template_dof; /**< 参考自由度. */
    BasisFunctionAdmin<double, DIM, DIM> basis_function; /**< 基函数 */
    std::vector<TemplateElement<double, DIM, DIM> > template_element; /**< 参考单元. */
    FEMSpace<double, DIM> fem_space; /**< 有限元空间. */
    EasyMesh mesh;		/**< 网格. */

    std::string mesh_file;	/**< 网格文件名. */
    double t;			/**< 物理时间. */
    double dt;			/**< 时间步长. */
    double gamma_x;		/**< 参数. */
    double gamma_y;		/**< 参数. */
    double beta;		/**< 参数. */
    double omega;		/**< 角速度. */
    FEMFunction<double, DIM> phi_re; /**< 数值解实部. */
    FEMFunction<double, DIM> phi_im; /**< 数值解虚部. */
    
public:
    /** 
     * 构造函数.
     * 
     * @param file 网格文件. 
     */
    RBEC(const std::string& file);
    /** 
     * 析构函数.
     * 
     * 
     * @return 标准返回.
     */
    virtual ~RBEC();
public:
    /** 
     * 主流程.
     * 
     */
    void run();
    /** 
     * 时间步长发展.
     * 
     */
    void stepForward();
    /** 
     * 初始化.
     * 
     */
    void initialize();
    /** 
     * 准备初值.
     * 
     */
    void initialValue();
    /** 
     * 计算能量.
     * 
     * @param phi_re 数值解实部.
     * @param phi_im 数值解虚部.
     * @param algebric_accuracy 积分精度.
     * 
     * @return 全局能量.
     */
    double energy(FEMFunction<double, DIM>& phi_re, FEMFunction<double, DIM>& phi_im,  int algebric_accuracy);    
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////

