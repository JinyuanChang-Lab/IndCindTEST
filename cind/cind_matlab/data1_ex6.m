function [ x,y,z ] = data1_ex6( n,p,q,m,dep)

    % 创建学生 t 分布实例
    x = trnd(2, [n,p]);
    y = trnd(2, [n,q]);
    z = trnd(2, [n, m]);
    w = trnd(2, [n,p]);
    v = trnd(2, [n,q]);
    % 获取所有可能的两两列组合的索引
    combinations = nchoosek(1:m, 2);
    num_combinations = size(combinations, 1);

    % 对每一对组合的对应元素进行相乘，并将结果存储在 x 矩阵中
    for i = 1:num_combinations
        col1 = combinations(i, 1);
        col2 = combinations(i, 2);
        x(:, i) = z(:, col1) .* z(:, col2);
        % x(:, i + num_combinations) = sin(z(:, col1)) + cos(z(:, col2)) + z(:, col1).^2 + z(:, col2).^2;
        y(:, i) = z(:, col1) + z(:, col2);
    end

    % 添加噪声
    if dep ~= 0
        for i = 1:dep
            eps = trnd(1, [n,1]);
            x(:, i) = x(:, i) + eps + 3 * eps.^3;
            y(:, i) = y(:, i) + eps + 3 * eps.^3;
        end
    end

end







