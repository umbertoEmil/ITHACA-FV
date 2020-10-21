// The TC file must have each TC on a column
// This function will mean the temperatures on the rows
Eigen::MatrixXd Tdanieli;
Tdanieli = cnpy::load(Tdanieli, TCfile, "colMajor");
Info << "Thermocouples temperature matrix shape = " << Tdanieli.rows() << " x "
     << Tdanieli.cols() << endl;
//std::cout << "Tdanieli = " << Tdanieli << std::endl;
Tdanieli = Tdanieli.array() + 273.15; // C to K
Tdanieli = ITHACAutilities::colWiseMean(Tdanieli);
//std::cout << "Tdanieli = " << Tdanieli << std::endl;
label TCnumber = thermocouplesCellID.size();
std::string message = "The file " + TCfile + " must contain " + std::to_string(
                          TCnumber) + " elements. It has " + std::to_string(Tdanieli.size());
M_Assert(Tdanieli.size() == TCnumber,
         message.c_str());
Tmeas = danieliTCtoMyNumeration(Tdanieli);
