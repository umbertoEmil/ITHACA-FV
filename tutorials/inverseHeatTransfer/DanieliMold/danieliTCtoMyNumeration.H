Eigen::MatrixXd Tout = Tdanieli;
label rows = 10;
label cols = 9;
M_Assert(Tdanieli.size() == rows * cols,
"Reprogram the function danieliTCtoMyNumeration with the right TC numbers.");
label i = 0;
for(int col = 0; col < cols; col++)
{
    for(int row = 0; row < rows; row++)
    {
        //Info << i << " = " << (rows - 1 - row) * cols + col << endl;
        Tout(i) = Tdanieli((rows - 1 - row) * cols + col);
        i++;
    }
}
return Tout;
