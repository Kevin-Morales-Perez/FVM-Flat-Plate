function [A_error_percentage] = matrix_err_percentage(A,A_correct)
%Evaluate the percentage of error of a matrix against a matrix of the same
%size with correct values 
if size(A)==size(A_correct)
    A_error_percentage=(abs(A-A_correct)./A_correct)*100;
else
    fprintf("Matrices are different size")

end

