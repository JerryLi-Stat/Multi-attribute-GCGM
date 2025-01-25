// [[Rcpp::depends(RcppEigen)]]
// #include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include <ctime>

using namespace std;

// Helper -----------------------------------------------------------------------------------

// Compute the sample covariance matrix of 'X'
Eigen::MatrixXd Covariance(const Eigen::MatrixXd &X)
{
    Eigen::MatrixXd X_center = X.array().rowwise() - X.array().colwise().mean();
    Eigen::MatrixXd cov = (X_center.transpose() * X_center) / (X_center.rows() - 1);
    return cov;
}

// [[Rcpp::export]]
Eigen::MatrixXd Soft_Threshold(const Eigen::MatrixXd &mat, double lambda)
{
    return mat.array().sign() * (mat.array().abs() - lambda).cwiseMax(0);
}

int num_unique(const Eigen::VectorXi &vec)
{
    std::unordered_set<int> unique_elements;

    for (int i = 0; i < vec.size(); i++)
    {
        unique_elements.insert(vec(i));
    }

    return unique_elements.size();
}

// 提取矩阵中指定行和列的元素
Eigen::MatrixXd extract_rows_cols(
    const Eigen::MatrixXd &mat,
    const vector<int> &rows, const vector<int> &cols)
{
    Eigen::MatrixXd sub_mat(rows.size(), cols.size());

    for (size_t i = 0; i < rows.size(); i++)
    {
        for (size_t j = 0; j < cols.size(); j++)
        {
            sub_mat(i, j) = mat(rows[i], cols[j]);
        }
    }

    return sub_mat;
}

// 修改矩阵中指定行和列的元素
void modify_matrix(Eigen::MatrixXd &mat, const std::vector<int> &rows, const std::vector<int> &cols, const Eigen::MatrixXd &new_values)
{
    for (size_t i = 0; i < rows.size(); ++i)
    {
        for (size_t j = 0; j < cols.size(); ++j)
        {
            mat(rows[i], cols[j]) = new_values(i, j);
        }
    }
}

Eigen::MatrixXi Find_Edge(
    const Eigen::MatrixXd &Omega,
    const Eigen::VectorXi &group_index)
{
    int n_groups = num_unique(group_index);
    Eigen::MatrixXi edge_fin(n_groups, n_groups);
    edge_fin.diagonal().setConstant(0);

    // the variables that contained in each group
    // assume group index starts from 1
    vector<vector<int>> attributes_group(n_groups);
    for (int i = 0; i < group_index.size(); i++)
    {
        attributes_group[group_index[i] - 1].push_back(i);
    }

    for (int i = 0; i < n_groups - 1; i++)
    {
        for (int j = i + 1; j <= n_groups - 1; j++)
        {
            double norm_val = extract_rows_cols(Omega, attributes_group[i], attributes_group[j]).norm();
            edge_fin(i, j) = (norm_val > 0) ? 1 : 0;
            edge_fin(j, i) = edge_fin(i, j);
        }
    }

    return edge_fin;
}

// Group Graphical Lasso ------------------------------------------------------------------------------------------

// Solution of equation (respect to A):
// A + \lambda * \sum_{i,j} A^{(ij)} / \|A^{(ij)}\|_F = B
// [[Rcpp::export]]
Eigen::MatrixXd Soft_threshold_Multi_attribute(
    const Eigen::MatrixXd &B, const Eigen::VectorXi &group_index,
    const double lambda, const bool diag_penalty = false)
{
    int p = B.cols();
    int n_groups = group_index.maxCoeff();
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(p, p);

    // the variables that contained in each group
    // assume group_index starts from 1
    vector<vector<int>> attributes_group(n_groups);
    for (int i = 0; i < group_index.size(); i++)
    {
        attributes_group[group_index[i] - 1].push_back(i);
    }

    // off-diagonal blocks
    for (size_t i = 0; i < n_groups - 1; i++)
    {
        for (size_t j = i + 1; j < n_groups; j++)
        {
            Eigen::MatrixXd B_ij = extract_rows_cols(B, attributes_group[i], attributes_group[j]);
            double c_ij = max(0.0, 1.0 - lambda / B_ij.norm());
            modify_matrix(A, attributes_group[i], attributes_group[j], B_ij * c_ij);
            modify_matrix(A, attributes_group[j], attributes_group[i], B_ij.transpose() * c_ij);
        }
    }

    // diagonal blocks
    for (size_t i = 0; i < n_groups; i++)
    {
        Eigen::MatrixXd B_ii = extract_rows_cols(B, attributes_group[i], attributes_group[i]);
        double c_ii = diag_penalty ? max(0.0, 1.0 - lambda / B_ii.norm()) : 1.0;
        modify_matrix(A, attributes_group[i], attributes_group[i], B_ii * c_ii);
    }

    return A;
}

Eigen::VectorXd Lambda_Select_GGL(
    const Eigen::MatrixXd &S, const Eigen::VectorXi &group_index,
    const int nLambda = 10, const double lambda_min_ratio = 0.1,
    const double rho_0 = 2.0)
{
    int n_groups = num_unique(group_index);

    // the variables that contained in each group
    // assume group_index starts from 1
    vector<vector<int>> attributes_group(n_groups);
    for (int i = 0; i < group_index.size(); i++)
    {
        attributes_group[group_index[i] - 1].push_back(i);
    }

    // Compute the smallest lambda that makes the first iteration returns a no-edge model
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_S(S);
    Eigen::MatrixXd V_S = eigen_S.eigenvectors();
    Eigen::VectorXd D_S = eigen_S.eigenvalues();
    Eigen::MatrixXd Omega_first_iter = V_S * (-D_S.array() + (D_S.array().square() + 4 * rho_0).sqrt()).matrix().asDiagonal() / (2 * rho_0) * V_S.transpose();

    double lambda_max = 0.0;
    for (size_t i = 0; i < n_groups - 1; i++)
    {
        for (size_t j = i + 1; j < n_groups; j++)
        {
            Eigen::MatrixXd B_ij = extract_rows_cols(Omega_first_iter, attributes_group[i], attributes_group[j]);
            lambda_max = max(B_ij.norm(), lambda_max);
        }
    }

    // The lambda_seq
    double lambda_min = lambda_max * lambda_min_ratio;
    Eigen::VectorXd lambdas = (Eigen::VectorXd::LinSpaced(nLambda, log(lambda_min), log(lambda_max))).array().exp();

    return lambdas.reverse();
}

// [[Rcpp::export]]
Eigen::MatrixXd GGL_ADMM(
    const Eigen::MatrixXd &X, const Eigen::VectorXi &group_index, double lambda,
    const Eigen::MatrixXd &W_init, const double rho_0 = 2.0, const double mu = 10.0,
    const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int iter_max = 1e3)
{
    int n = X.rows(), p = X.cols();
    int n_groups = num_unique(group_index);

    // the variables that contained in each group
    // assume group_index starts from 1
    vector<vector<int>> attributes_group(n_groups);
    for (int i = 0; i < group_index.size(); i++)
    {
        attributes_group[group_index[i] - 1].push_back(i);
    }

    Eigen::MatrixXd S = Covariance(X);
    Eigen::MatrixXd U_old = Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd Omega_old = S.diagonal().cwiseInverse().asDiagonal();
    Eigen::MatrixXd W_old = Eigen::MatrixXd::Zero(p, p);
    if (W_init.rows() == p && W_init.cols() == p)
    {
        W_old = W_init;
    }
    Eigen::MatrixXd W_new = Eigen::MatrixXd::Zero(p, p);

    double rho = rho_0;
    bool converged = false;
    int iter = 0;

    while (!converged && (iter <= iter_max))
    {
        // Update Omega
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_temp(S - rho * (W_old - U_old));
        Eigen::MatrixXd V_temp = eigen_temp.eigenvectors();
        Eigen::VectorXd D_temp = eigen_temp.eigenvalues();
        Eigen::MatrixXd D_tilde = ((-D_temp.array() + (D_temp.array().square() + 4 * rho).sqrt()) / (2 * rho)).matrix().asDiagonal();
        Eigen::MatrixXd Omega_new = V_temp * D_tilde * V_temp.transpose();

        // Update W
        W_new = Soft_threshold_Multi_attribute(Omega_new + U_old, group_index, lambda / rho);

        // Update U
        Eigen::MatrixXd U_new = U_old + Omega_new - W_new;

        // Stopping criterion
        double d_p = (Omega_new - W_new).norm();
        double d_d = rho * (Omega_new - Omega_old).norm();
        double tol_pri = p * tol_abs + tol_rel * max(Omega_new.norm(), W_new.norm());
        double tol_dual = p * tol_abs + tol_rel * U_new.norm() * rho; // it should be ||U||_F * rho here

        if (d_p <= tol_pri && d_d <= tol_dual)
        {
            converged = true;
        }

        if (d_p > mu * d_d)
        {
            rho *= 2;
            U_new /= 2;
        }

        if (d_d > mu * d_p)
        {
            rho /= 2;
            U_new *= 2;
        }

        W_old = W_new;
        Omega_old = Omega_new;
        U_old = U_new;
        ++iter;
    }

    return W_new;
}

// [[Rcpp::export]]
Rcpp::List GGL_BIC(
    const Eigen::MatrixXd &X, const Eigen::VectorXi &group_index,
    const int nLambda = 10, const double lambda_min_ratio = 0.1,
    const double rho_0 = 2.0, const double mu = 10.0,
    const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int iter_max = 1e3)
{
    clock_t start = clock();
    int n = X.rows(), p = X.cols();
    int n_groups = num_unique(group_index);
    Eigen::MatrixXd S = Covariance(X);

    // lambda_seq
    Eigen::VectorXd lambda_seq = Lambda_Select_GGL(S, group_index, nLambda, lambda_min_ratio, rho_0);

    // Select lambda
    vector<Eigen::MatrixXd> Omega_list(nLambda);
    Eigen::VectorXd BICs(nLambda);

    for (size_t l = 0; l < nLambda; ++l)
    {
        double lambda_l = lambda_seq(l);
        Eigen::MatrixXd W_init = (l == 0 ? Eigen::MatrixXd::Zero(p, p) : Omega_list[l - 1]);
        Eigen::MatrixXd Omega_l = GGL_ADMM(X, group_index, lambda_l, W_init, rho_0, mu, tol_abs, tol_rel, iter_max);
        double BIC_l = (S * Omega_l).trace() - log(Omega_l.determinant()) + log(n) * (Omega_l.array() != 0).count() / (2 * n);
        Omega_list[l] = Omega_l;
        BICs(l) = BIC_l;
    }

    Eigen::Index min_lambda_idx;
    double minBIC = BICs.minCoeff(&min_lambda_idx);
    double lambda_final = lambda_seq(min_lambda_idx);
    Eigen::MatrixXd Omega_final = Omega_list[min_lambda_idx];

    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;

    auto Edge_final = Find_Edge(Omega_final, group_index);

    return Rcpp::List::create(
        Rcpp::Named("Omega") = Omega_final, Rcpp::Named("Edge") = Edge_final,
        Rcpp::Named("lambda_seq") = lambda_seq, Rcpp::Named("lambda_final") = lambda_final,
        Rcpp::Named("time") = duration);
}

// Sparse-Group Graphical Lasso -----------------------------------------------------------------------------------

// [[Rcpp::export]]
Eigen::VectorXd Lambda_Select_SGL(
    const Eigen::MatrixXd &X, const Eigen::VectorXi &group_index, const Eigen::MatrixXd &S,
    int nLambda, double alpha_0 = 0.1, double rho_0 = 2.0,
    double lambda_min_ratio = 0.1, double tau = 0.002)
{
    int n_groups = num_unique(group_index);
    alpha_0 = min(max(alpha_0, 0.1), 0.9);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_S(S);
    Eigen::MatrixXd V_S = eigen_S.eigenvectors();
    Eigen::VectorXd D_S = eigen_S.eigenvalues();
    Eigen::MatrixXd Omega_first_iter = V_S * (-D_S.array() + (D_S.array().square() + 4 * rho_0).sqrt()).matrix().asDiagonal() / (2 * rho_0) * V_S.transpose();
    Omega_first_iter.diagonal().setZero();
    double lambda_extrem = Omega_first_iter.cwiseAbs().maxCoeff() * rho_0 / alpha_0;
    double lambda_max = lambda_extrem / 2;
    double lambda_min = lambda_max * lambda_min_ratio;
    Eigen::VectorXd lambdas = (Eigen::VectorXd::LinSpaced(nLambda, log(lambda_min), log(lambda_max))).array().exp();

    return lambdas.reverse();
}

// [[Rcpp::export]]
Eigen::MatrixXd SGL_ADMM(
    const Eigen::MatrixXd &X, const Eigen::VectorXi &group_index, double lambda, double alpha,
    const Eigen::MatrixXd &W_init, const int cov_method = 1,
    const double rho_0 = 2.0, const double mu = 10.0,
    const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int iter_max = 1e3)
{
    int n = X.rows(), p = X.cols();
    int n_groups = num_unique(group_index);

    // the variables that contained in each group
    // assume group_index starts from 1
    vector<vector<int>> attributes_group(n_groups);
    for (int i = 0; i < group_index.size(); i++)
    {
        attributes_group[group_index[i] - 1].push_back(i);
    }

    Eigen::MatrixXd X_centered = X.rowwise() - X.colwise().mean();
    Eigen::MatrixXd S(p, p);
    if (cov_method == 1)
    {
        S = X_centered.transpose() * X_centered / n;
    }
    else if (cov_method == 2)
    {
        Rcpp::Function cor_fk("cor.fk", Rcpp::Environment::namespace_env("pcaPP"));
        Rcpp::NumericMatrix cor_matrix = cor_fk(Rcpp::wrap(X));
        S = (Rcpp::as<Eigen::MatrixXd>(cor_matrix).array() * M_PI / 2).sin();
    }

    Eigen::MatrixXd U_old = Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd Omega_old = S.diagonal().cwiseInverse().asDiagonal();
    Eigen::MatrixXd W_old = Eigen::MatrixXd::Zero(p, p);
    if (W_init.rows() == p && W_init.cols() == p)
    {
        W_old = W_init;
    }
    Eigen::MatrixXd W_new = Eigen::MatrixXd::Zero(p, p);

    double rho = rho_0;
    bool converged = false;
    int iter = 0;

    while (!(converged) && (iter <= iter_max))
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_temp(S - rho * (W_old - U_old));
        Eigen::MatrixXd V_temp = eigen_temp.eigenvectors();
        Eigen::VectorXd D_temp = eigen_temp.eigenvalues();
        Eigen::MatrixXd D_tilde = ((-D_temp.array() + (D_temp.array().square() + 4 * rho).sqrt()) / (2 * rho)).matrix().asDiagonal();
        Eigen::MatrixXd Omega_new = V_temp * D_tilde * V_temp.transpose();

        W_new.setZero();
        Eigen::MatrixXd A_temp = Omega_new + U_old;

        // off-diagonal subblocks of W
        for (size_t i = 0; i < n_groups - 1; i++)
        {
            for (size_t j = i + 1; j < n_groups; j++)
            {
                Eigen::MatrixXd B_temp = Soft_Threshold(extract_rows_cols(A_temp, attributes_group[i], attributes_group[j]), alpha * lambda / rho);

                modify_matrix(W_new, attributes_group[i], attributes_group[j], B_temp * max(0.0, 1 - lambda * (1 - alpha) / (rho * B_temp.norm())));
                modify_matrix(W_new, attributes_group[j], attributes_group[i], B_temp.transpose() * max(0.0, 1 - lambda * (1 - alpha) / (rho * B_temp.norm())));
            }
        }
        // W_new += W_new.transpose();

        // diagonal subblocks of W
        for (size_t i = 0; i < n_groups; i++)
        {
            Eigen::MatrixXd A_ii = extract_rows_cols(A_temp, attributes_group[i], attributes_group[i]);
            Eigen::MatrixXd Diag_of_Aii = A_ii.diagonal().asDiagonal();
            modify_matrix(W_new, attributes_group[i], attributes_group[i], Soft_Threshold(A_ii - Diag_of_Aii, alpha * lambda / rho) + Diag_of_Aii);
        }

        Eigen::MatrixXd U_new = U_old + Omega_new - W_new;
        double d_p = (Omega_new - W_new).norm();
        double d_d = rho * (Omega_new - Omega_old).norm();
        double tol_pri = p * tol_abs + tol_rel * max(Omega_new.norm(), W_new.norm());
        double tol_dual = p * tol_abs + tol_rel * U_new.norm() * rho; // it should be ||U||_F * rho here

        if (d_p <= tol_pri && d_d <= tol_dual)
        {
            converged = true;
        }

        if (d_p > mu * d_d)
        {
            rho *= 2;
            U_new /= 2;
        }

        if (d_d > mu * d_p)
        {
            rho /= 2;
            U_new *= 2;
        }

        W_old = W_new;
        Omega_old = Omega_new;
        U_old = U_new;
        ++iter;
    }

    return W_new;
}

// [[Rcpp::export]]
Rcpp::List SGL_BIC(
    const Eigen::MatrixXd &X, const Eigen::VectorXi &group_index,
    const int nLambda = 10, const double lambda_min_ratio = 0.1,
    const double alpha_0 = 0.1, const int nalpha = 9, const double alpha_min = 0.1, const double alpha_max = 0.9,
    const int cov_method = 1, const bool use_corr = true,
    const double rho_0 = 2.0, const double mu = 10.0,
    const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int iter_max = 1e3)
{
    clock_t start = clock();
    int n = X.rows(), p = X.cols();
    int n_groups = num_unique(group_index);

    // Compute Sigma_hat
    Eigen::MatrixXd X_centered = X.rowwise() - X.colwise().mean();
    Eigen::VectorXd X_sd = ((X_centered.array().square().colwise().sum()) / (n - 1)).sqrt();
    Eigen::MatrixXd S(p, p);

    if (cov_method == 1)
    {
        S = X_centered.transpose() * X_centered / n;
        if (use_corr)
        {
            for (int i = 0; i < p; ++i)
            {
                X_centered.col(i) /= X_sd(i);
            }

            S = X_sd.cwiseInverse().asDiagonal() * S * X_sd.cwiseInverse().asDiagonal();
        }
    }
    else if (cov_method == 2)
    {
        Rcpp::Function cor_fk("cor.fk", Rcpp::Environment::namespace_env("pcaPP"));
        Rcpp::NumericMatrix cor_matrix = cor_fk(Rcpp::wrap(X));
        S = (Rcpp::as<Eigen::MatrixXd>(cor_matrix).array() * M_PI / 2).sin();
    }

    // Compute lambdas
    Eigen::VectorXd lambda_seq = Lambda_Select_SGL(X_centered, group_index, S, nLambda, alpha_0, rho_0, lambda_min_ratio);
    Eigen::VectorXd alpha_seq = Eigen::VectorXd::LinSpaced(nalpha, alpha_min, alpha_max);

    // Select lambda with alpha = alpha_0
    vector<Eigen::MatrixXd> Omega_list(nLambda);
    Eigen::VectorXd BIC_lambda(nLambda);

    for (size_t l = 0; l < nLambda; ++l)
    {
        double lambda_l = lambda_seq(l);
        Eigen::MatrixXd W_init = (l == 0 ? Eigen::MatrixXd::Zero(p, p) : Omega_list[l - 1]);
        Eigen::MatrixXd Omega_l = SGL_ADMM(X_centered, group_index, lambda_l, alpha_0, W_init, cov_method,
                                           rho_0, mu, tol_abs, tol_rel, iter_max);
        double BIC_l = (S * Omega_l).trace() - log(Omega_l.determinant()) + log(n) * (Omega_l.array() != 0).count() / (2 * n);
        Omega_list.push_back(Omega_l);
        BIC_lambda(l) = BIC_l;
    }

    Eigen::Index min_lambda_idx;
    double minBIC = BIC_lambda.minCoeff(&min_lambda_idx);
    double lambda_final = lambda_seq(min_lambda_idx);

    // Select alpha with lambda = lambda_final
    vector<Eigen::MatrixXd>(nalpha).swap(Omega_list);
    Eigen::VectorXd BIC_alpha(nalpha);

    for (size_t l = 0; l < nalpha; ++l)
    {
        double alpha_l = alpha_seq(l);
        Eigen::MatrixXd W_init = (l == 0 ? Eigen::MatrixXd::Zero(p, p) : Omega_list[l - 1]);
        Eigen::MatrixXd Omega_l = SGL_ADMM(X_centered, group_index, lambda_final, alpha_l, W_init, cov_method,
                                           rho_0, mu, tol_abs, tol_rel, iter_max);
        double BIC_l = (S * Omega_l).trace() - log(Omega_l.determinant()) + log(n) * (Omega_l.array() != 0).count() / (2 * n);
        Omega_list[l] = Omega_l;
        BIC_alpha(l) = BIC_l;
    }

    Eigen::Index min_alpha_idx;
    minBIC = BIC_alpha.minCoeff(&min_alpha_idx);
    double alpha_final = alpha_seq(min_alpha_idx);
    Eigen::MatrixXd Omega_final = Omega_list[min_alpha_idx];

    // Adjust the final result
    if (cov_method == 2 || use_corr)
    {
        Omega_final = X_sd.cwiseInverse().asDiagonal() * Omega_final * X_sd.cwiseInverse().asDiagonal();
    }

    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;

    auto Edge_final = Find_Edge(Omega_final, group_index);

    return Rcpp::List::create(
        Rcpp::Named("Omega") = Omega_final, Rcpp::Named("Edge") = Edge_final,
        Rcpp::Named("lambda_seq") = lambda_seq, Rcpp::Named("lambda_final") = lambda_final,
        Rcpp::Named("alpha_seq") = alpha_seq, Rcpp::Named("alpha_final") = alpha_final,
        Rcpp::Named("time") = duration);
}
