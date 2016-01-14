/**
 * @file   BEC.h
 * @author Feng Yue <fengyue@zju.edu.cn>
 * @date   Wed Dec 23 09:58:53 2015
 * 
 * @brief  BEC方程的类型声明.
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
    double omega0;		/**< 参数. */
    double delta;		/**< 参数. */
    double r0;		        /**< 参数. */
public:
    /** 
     * 初始构造函数, 传递参数.
     * 
     * @param _gamma_x 参数. 
     * @param _gamma_y 参数.
     * @param _omega0  参数.
     * @param _delta   参数.
     * @param _r0      参数.
     * 
     * @return 标准返回.
     */
    Potential(double _gamma_x, 
	      double _gamma_y,
	      double _omega0,
	      double _delta,
              double _r0) :
	gamma_x(_gamma_x), 
	gamma_y(_gamma_y),
	omega0(_omega0),
        delta(_delta),
        r0(_r0)
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
 *  初值函数类型.
 * 
 */
class Initial : public Function<double>
{
private:
    double gamma_x;		/**< 参数. */
    double gamma_y;		/**< 参数. */
public:
    /** 
     * 初始构造函数.
     * 
     * @param _gamma_x 参数.
     * @param _gamma_y 参数.
     * 
     * @return 标准返回.
     */
    Initial(double _gamma_x, 
	    double _gamma_y) :
	gamma_x(_gamma_x), 
	gamma_y(_gamma_y)
	{};
    /** 
     * 析构函数.
     * 
     * 
     * @return 标准返回. 
     */
    ~Initial() {};
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
 * BEC方程类.
 * 
 */
class BEC
{
public:
    /**
     * 系数矩阵类.
     * 
     */
    class Matrix : public StiffMatrix<DIM, double>
    {
    private:
	double dt;		/**< 时间步长. */
	double beta;		/**< 参数. */
	double gamma_x;		/**< 参数. */
	double gamma_y;		/**< 参数. */
	double omega0;		/**< 参数. */
	double delta;		/**< 参数. */
	double r0;		/**< 参数. */
	FEMFunction<double, DIM> *phi; /**< 数值解. */
    public:
	/** 
	 * 构造函数.
	 * 
	 * @param sp 有限元空间.
	 * @param _dt 时间步长.
	 * @param _beta 参数.
	 * @param _gamma_x 参数.
	 * @param _gamma_y 参数.
	 * @param _omega0  参数.
	 * @param _delta   参数.
	 * @param _r0      参数.
	 * @param _phi 数值解, 实部.
	 * 
	 * @return 标准返回.
	 */
	Matrix(FEMSpace<double, DIM>& sp, 
		 const double& _dt, 
		 const double& _beta,
		 const double& _gamma_x,
		 const double& _gamma_y,
                 const double& _omega0,
                 const double& _delta,
                 const double& _r0,
		 FEMFunction<double, DIM>& _phi) :
	    dt(_dt), 
	    beta(_beta),	
	    gamma_x(_gamma_x),
	    gamma_y(_gamma_y),
            omega0(_omega0),
            delta(_delta),
            r0(_r0),
	    StiffMatrix<DIM, double>(sp) 
	    {
		phi = &_phi;
	    };
	/** 
	 * 析构函数.
	 * 
	 * 
	 * @return 标准返回. 
	 */
	virtual ~Matrix() {};
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
    double omega0;		/**< 参数. */
    double delta; 		/**< 参数. */
    double r0;  		/**< 参数. */
    double beta;		/**< 参数. */
    FEMFunction<double, DIM> phi; /**< 数值解实部. */
    
public:
    /** 
     * 构造函数.
     * 
     * @param file 网格文件. 
     */
    BEC(const std::string& file);
    /** 
     * 析构函数.
     * 
     * 
     * @return 标准返回.
     */
    virtual ~BEC();
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
     * @param phi 数值解实部.
     * @param algebric_accuracy 积分精度.
     * 
     * @return 全局能量.
     */
    double energy(FEMFunction<double, DIM>& phi,  int algebric_accuracy);    
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////

