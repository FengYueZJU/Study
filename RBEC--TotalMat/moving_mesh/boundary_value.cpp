#include "RBEC.h"
#define DIM 2

void RBEC::boundaryValue(Vector<double> &x, Vector<double> &rhs, SparseMatrix<double> &matrix)
{
	/// 空间自由度.
	unsigned int n_dof = fem_space.n_dof();
	unsigned int n_total_dof = 2 * n_dof;
	const std::size_t * rowstart = sp_RBEC.get_rowstart_indices();
	const unsigned int * colnum = sp_RBEC.get_column_numbers();
//	std::cout << "n_dof: " << n_dof << std::endl;
//	std::cout << "n_A: " << sp_RBEC.n_rows() << ", m_A: " << sp_RBEC.n_cols() << std::endl;

	/// 遍历全部维度的节点.
	for (unsigned int i = 0; i < n_total_dof; ++i)
	{
		/// 边界标志.
		int bm = -1;
		/// 判断一下是实部还是虚部，分别读取标志.
		if (i < n_dof)
			bm = fem_space.dofInfo(i).boundary_mark;
		else
			bm = fem_space.dofInfo(i - n_dof).boundary_mark;

		if (bm == 0)
			continue;

		/// 零边界条件.
		if (bm == 1)
			x(i) = 0.0;

		/// 右端项这样改, 如果该行和列其余元素均为零, 则在迭代中确保该数值解和边界一致.

		if (bm == 1)
		{
			rhs(i) = matrix.diag_element(i) * x(i);
			/// 遍历 i 行.
			for (unsigned int j = rowstart[i] + 1;j < rowstart[i + 1]; ++j)
			{
				/// 第 j 个元素消成零(不是第 j 列!). 注意避开了对角元.
				matrix.global_entry(j) -= matrix.global_entry(j);
				/// 第 j 个元素是第 k 列.
				unsigned int k = colnum[j];
				/// 看看第 k 行的 i 列是否easymesh 为零元.
				const unsigned int *p = std::find(&colnum[rowstart[k] + 1],&colnum[rowstart[k + 1]],i);
				/// 如果是非零元. 则需要将这一项移动到右端项. 因为第 i 个未知量已知.
				if (p != &colnum[rowstart[k + 1]])
				{
					/// 计算 k 行 i 列的存储位置.
					unsigned int l = p - &colnum[rowstart[0]];
					/// 移到右端项. 等价于 r(k) = r(k) - x(i) * A(k, i).
					rhs(k) -= matrix.global_entry(l)* x(i);
					/// 移完此项自然是零.
					matrix.global_entry(l) -= matrix.global_entry(l);
				}
			}
		}
	}
//	std::cout << "Boundary values for RBEC OK!" << std::endl;
};

#undef DIM
