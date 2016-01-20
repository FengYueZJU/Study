#include "RBEC.h"
#define DIM 2

void RBEC::buildMatrixStruct()
{
    /// 构建系数矩阵模板.
    /// 空间自由度.
    int n_dof = fem_space.n_dof();
    /// 总自由度.
    int n_total_dof = 2 * n_dof;
    /// 准备统计系数矩阵的每一行有多少个非零元.
    std::vector<unsigned int> n_non_zero_per_row(n_total_dof);
    /// 设置各块大小.
    int row_re = n_dof;	/**< (0, 0) */
    int col_re = n_dof;
    std::vector<unsigned int> n_non_zero_per_row_rere(n_dof);
    int col_im = n_dof;       /**< (0, 1) */
    std::vector<unsigned int> n_non_zero_per_row_reim(n_dof);
    int row_im = n_dof;	/**< (1, 0) */
    std::vector<unsigned int> n_non_zero_per_row_imre(n_dof);
    /**< (1, 1) */
    std::vector<unsigned int> n_non_zero_per_row_imim(n_dof);

    /// 准备一个遍历全部单元的迭代器.
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space.endElement();

    /// 第一次循环遍历全部单元, 只是为了统计每一行的非零元个数.
    for (; the_element != end_element; ++the_element)
    {
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
	/// j 是检测函数指标. 行指标. 在整个计算中, 默认行是检测函数而
	/// 列是基函数.
	for (int j = 0; j < n_element_dof; ++j)
	{
	    /// k 是基函数指标. 列指标.
	    for (int k = 0; k < n_element_dof; ++k)
	    {
		/// element_dof[j] 行两个.
		n_non_zero_per_row[element_dof[j]] += 2;
		/// 一个在 (0, 0) 块.
		n_non_zero_per_row_rere[element_dof[j]]++;
		/// 一个在 (0, 1) 块.
		n_non_zero_per_row_reim[element_dof[j]]++;
		/// n_dof_v + element_dof[j] 行两个.
		n_non_zero_per_row[element_dof[j] + n_dof] += 2;
		/// 这个在 (1, 0) 块.
		n_non_zero_per_row_imre[element_dof[j]]++;
		/// 这个在 (1, 1) 块.
		n_non_zero_per_row_imim[element_dof[j]]++;
	    }
	}
    }

    /// 指定矩阵模板带宽.
    sp_RBEC.reinit(n_total_dof, n_total_dof, n_non_zero_per_row, true);
    /// 应用各块带宽.
    sp_rere.reinit(row_re, col_re, n_non_zero_per_row_rere, true);	/**< (0, 0) */
    sp_reim.reinit(row_re, col_im, n_non_zero_per_row_reim, false);	/**< (0, 1) */
    sp_imre.reinit(row_im, col_re, n_non_zero_per_row_imre, false);	/**< (1, 0) */
    sp_imim.reinit(row_re, col_im, n_non_zero_per_row_imim, true);	/**< (1, 1) */

    /// 第二次遍历, 指定每个非零元的坐标.
    for (the_element = fem_space.beginElement();
	 the_element != end_element; ++the_element)
    {
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
	/// j 是检测函数指标. 行指标.
	for (int j = 0; j < n_element_dof; ++j)
	{
	    for (int k = 0; k < n_element_dof; ++k)
	    {
		/// (0, 0)
		sp_RBEC.add(element_dof[j], element_dof[k]);
		sp_rere.add(element_dof[j], element_dof[k]);
		/// (0, 1)
		sp_RBEC.add(element_dof[j], element_dof[k] + n_dof);
		sp_reim.add(element_dof[j], element_dof[k]);
		/// (1, 0)
		sp_RBEC.add(element_dof[j] + n_dof, element_dof[k]);
		sp_imre.add(element_dof[j], element_dof[k]);
		/// (1, 1)
		sp_RBEC.add(element_dof[j] + n_dof, element_dof[k] + n_dof);
		sp_imim.add(element_dof[j], element_dof[k]);
	    }
	}
    }

    /// 矩阵模板压缩. 创建矩阵.
    sp_RBEC.compress();
    /// 压缩各块模板.
    sp_rere.compress();	/**< (0, 0) */
    sp_reim.compress();	/**< (0, 1) */
    sp_imre.compress();	/**< (1, 0) */
    sp_imim.compress();	/**< (1, 1) */

    index_rere.resize(sp_rere.n_nonzero_elements());
    index_imre.resize(sp_imre.n_nonzero_elements());
    index_reim.resize(sp_reim.n_nonzero_elements());
    index_imim.resize(sp_imim.n_nonzero_elements());

    std::vector<int>::iterator index_rere_iterator = index_rere.begin();
    std::vector<int>::iterator index_reim_iterator = index_reim.begin();
    std::vector<int>::iterator index_imre_iterator = index_imre.begin();
    std::vector<int>::iterator index_imim_iterator = index_imim.begin();

    /// 准备整体矩阵索引.
    const std::size_t * rowstart = sp_RBEC.get_rowstart_indices();
    const unsigned int * column = sp_RBEC.get_column_numbers();
    /// 准备各块列索引.
    const unsigned int * column_rere = sp_rere.get_column_numbers();
    const unsigned int * column_reim = sp_reim.get_column_numbers();
    const unsigned int * column_imre = sp_imre.get_column_numbers();
    const unsigned int * column_imim = sp_imim.get_column_numbers();

    for (int i = 0; i < n_total_dof; ++i)
    {
	int l_i1 = row_re;
	int l_j1 = col_re;
	for (int j = rowstart[i]; j < rowstart[i + 1]; ++j)
	{
	    if (i >= 0 && i < l_i1)
	    {
		if (column[j] >= 0 && column[j] < l_j1)
		{
		    if (column[j] !=  *column_rere)
		    {
			std::cout << "Matrix structure error, (0, 0)!" << std::endl;
			exit(-1);
		    }
		    else
		    {
			*index_rere_iterator = j;
			index_rere_iterator++;
			column_rere++;
		    }
		}
		else if (column[j] >= l_j1)
		{
		    if (column[j] - l_j1 !=  *column_reim)
		    {
			std::cout << "Matrix structure error, (0, 1)!" << std::endl;
			exit(-1);
		    }
		    else
		    {
			*index_reim_iterator = j;
			index_reim_iterator++;
			column_reim++;
		    }
		}
	    }
	    else if (i >= l_i1)
	    {
		if (column[j] >= 0 && column[j] < l_j1)
		{
		    if (column[j] !=  *column_imre)
		    {
			std::cout << "Matrix structure error, (1, 0)!!" << std::endl;
			exit(-1);
		    }
		    else
		    {
			*index_imre_iterator = j;
			index_imre_iterator++;
			column_imre++;
		    }
		}
		else if (column[j] >= l_j1)
		{
		    if (column[j] - l_j1 !=  *column_imim)
		    {
			std::cout << "Matrix structure error, (1, 1)!" << std::endl;
			exit(-1);
		    }
		    else
		    {
			*index_imim_iterator = j;
			index_imim_iterator++;
			column_imim++;
		    }
		}
	    }
	}
    }
    std::cout << "Matrix struct partion OK!" << std::endl;
};

#undef DIM
