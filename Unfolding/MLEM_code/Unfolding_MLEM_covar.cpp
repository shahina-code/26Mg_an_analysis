/************************************************************************
 *
 *  Filename: Unfolding.cpp
 *
 *  Description:
 *    Maximum-Likelihood Expectation Maximization (MLEM)
 *
 *  Compile:
 *	  To compile use "g++ Unfolding.cpp -I. -o Program.exe"
 *
 *  Inputs:
 *    spectrum:   spectrum to be unfolded
 *    iterations: number of interations*
 *    binWidth: width of bins in spectrum in MeVee
 *    options:
 *          -c  terminate on convergence (default)
 *          -i  terminate when iterations complete
 *          -p  print result to output file
 *
 *
 *  Author(s):
 *     Haonan Zhu
 *     Michael T. Febbraro
 *
 *
 *  Creation Date: 2/8/2014
 *  Last modified: 2/8/2014
 *
 * -----------------------------------------------------
 *  Nuclear Reaction Group
 *  University of Michigan, Ann Arbor, MI, USA
 *  (c) All Rights Reserved.
 *
 */

#include <iostream>
#include <sstream>
#include <string>
#include "/Users/shahinashahina/Documents/PhD/Experiments/25Mg(a,n)/CAMPAIGN_2021/25Mg_codes_Shahina/eigen-3.3.7/Eigen/Eigen"
#include <fstream>
#include <queue>

using namespace Eigen;
using namespace std;

/** ----------------------------------------------------
* Variable declairation
*   ----------------------------------------------------
*/

  bool run = 1;
  double val;
  char filename[250], prompt;
  int I, J, k, iterations, eBins, lBins,sLength, NPS, N;
  int Nconstrains, Nmatrix[100], matrixIndex;
  double norm_factor, threshold, criteria, max_diff;
              //I: Number of inputs, J:Number of Outputs
  float lwidth, binWidth;
  MatrixXd R;  //Response Matrix
  MatrixXd V;  // Diagonized input
  MatrixXd y;  //Input
  MatrixXd Z; //Inverse of Response Matrix
  VectorXd Rs; //Vector of sums
  VectorXd x; //Output
  VectorXd diff; //difference between S(k+1) - S(k)
  VectorXd r; //Residuals
  VectorXd residual; //Residuals
  VectorXd dx; // Uncertainty
  VectorXd ax; // Average of x
  VectorXd t;  // temporary vector
  VectorXd s;  // Scaling vector
  VectorXd ToF;  // ToF vector
  VectorXd iter;
  VectorXd chi2;
  VectorXd IndicatorFunc;
  VectorXd IndicatorMaxDiff;
  VectorXd PeakSum, PeakSum2;
  MatrixXd sum;
  string line;
  clock_t Stop, Start;
  stringstream ss (stringstream::in | stringstream::out);
  ofstream out;

namespace Eigen

// Define a new fstream input method for matrix
{
    fstream& operator>>(fstream &fin, MatrixXd &mat) {
            int h=0, w,size=0;
            queue<double> In;
            getline(fin,line);
            ss.clear();
            ss<<line;
            for(size=0;ss>>val;size++);
            ss.clear();
            ss<<line;
            if(size == 3){
                 ss>>sLength;     cout << " - Number of bins: " << sLength << endl;
                 ss>>binWidth;    cout << " - Bin width: " << binWidth << " MeVee" << endl;
                 ss>>iterations;  cout << " - Iterations: " << iterations << endl;
            }
            else if(size == 4){
                  ss>> lBins; cout <<" - Light bins:" << lBins << endl;
                  ss >> eBins; cout <<" - Energy bins: " << eBins << endl;
                  ss>> NPS; cout << " - NPS: " << NPS << endl;
                  ss>> lwidth; cout << " - Bin width: " << lwidth << " MeVee" << endl;
            }
            else{
              cout<<"Wrong File Formatting ";
              exit(1);
            }
            while(getline(fin, line)){
                if(line.empty()) break;
                ss.clear();
                ss << line;
                while(ss >> val)
                {
			        //if(val<=threshold) val=0;
                    In.push(val);
                }
                if(h==0) w = In.size();
                h++;
            }
            mat = MatrixXd(h, w);
            int i=0,j=0;
            while(!In.empty()) {
                mat(i,j) = In.front();
                j++;
                if(j==w) {
                    j=0, i++;
                }
                In.pop();
            }
            return fin;

        }
}

MatrixXd Update(MatrixXd x,MatrixXd zero,int I)
//Replace all zeros terms in the matrix with 1
{
     MatrixXd check(I,1);
     stringstream ss (stringstream::in | stringstream::out);
     ss.clear();
     ss<<(x.array()==zero.array());
     queue<double> In;
     double val;
     while(ss >> val)
     {
       	In.push(val);
     }

     int i=0;
     while(i<I) {
         check(i,0) = In.front();
         i++;
         In.pop();
     }
     return x+check;
}

int main(int argc, char** argv)
{
    Start = clock();


	cout << "--------------------------------------" << endl;
	cout << "   MLEM spectrum unfolding            " << endl;
  cout << "   - Note: covariance matrix added    " << endl;
	cout << "--------------------------------------" << endl;
	cout << endl;

	//cout << "Spectrum File: ";
	//cin >> filename;

	//cout << "Threshold: ";
	//cin >> threshold;
  threshold = 20;

  criteria = 1.0E-6;
  cout << "Stopping criteria : " << criteria << endl;

	while (run)
	{

    /** ----------------------------------------------------
    * Load response matrix
    *   ----------------------------------------------------
    */
    cout << "\n" << "Loading response matrix..." << endl;

    //fstream fin("matrix.txt",fstream::in);
    fstream  fin("/Users/shahinashahina/Documents/PhD/Experiments/25Mg(a,n)/CAMPAIGN_2021/25Mg_codes_Shahina/matrix.txt",fstream::in);
    fin>>Z;
    cout << "Response Matrix loaded!\n" << endl;


    /** ----------------------------------------------------
    * Load spectrum
    *   ----------------------------------------------------
    */
    cout << "\n" << "Loading input spectrum..." << endl;
    //fstream fan(filename,fstream::in);
    fstream fan("spectrum.spe" ,fstream::in);
    fan>>y;
    cout << "Spectrum loaded!\n" << endl;

    I=(sLength>=lBins)*lBins+(sLength<lBins)*sLength;
    J=eBins;

    /** ----------------------------------------------------
    * Load ToF spectrum
    *   ----------------------------------------------------
    */
    ToF=VectorXd::Zero(J+1);
    cout << "\n" << "Loading tof spectrum..." << endl;
    fstream fann("ToF.spe",fstream::in);
    int jj = 0; double sum_tof = 0;
    while(getline(fann, line) && jj < J){
      ToF(jj) = atof(line.c_str());
      sum_tof += ToF(jj);
      jj++;
    }

    for (int jjj = 0; jjj < jj; jjj++)
    {
      ToF(jjj) = (ToF(jjj)/sum_tof)*y.sum();
    }

    cout << "Spectrum loaded!\n" << endl;

    // Override number of Iterations
    //iterations = atoi(argv[1]);
    cout << "interations : " << iterations << endl;

	/** ----------------------------------------------------
    * Select operation mode
    *   ----------------------------------------------------
    */
	//cout << "Constrain matrix? (y/n):";
	//cin >> prompt;
  prompt = 'n';
	s=VectorXd::Zero(J);
	if (prompt == 'y' || prompt == 'Y')
	{
		cout << "Number of contrains: ";
		cin >> Nconstrains;

		for (int i = 0; i < Nconstrains; i++)
		{
			cout << "Value: ";
			cin >> matrixIndex;
			s(matrixIndex) = 1;
		}

	}



    /** ----------------------------------------------------
    * normalize, & rebin response matrix & input Spectrum
    * initilize important operational matrices/vectors
    *   ----------------------------------------------------
    */

	for (int i = 0; i < y.size(); i++)
	{
		if(i < threshold) {y(i) = 0;}
	}

  x=VectorXd::Ones(J);
	dx=VectorXd::Zero(J);
	residual=VectorXd::Zero(iterations);
  iter = VectorXd::Zero(iterations);
  chi2 = VectorXd::Zero(iterations);
  IndicatorFunc = VectorXd::Zero(iterations);
  IndicatorMaxDiff = VectorXd::Zero(iterations);
  PeakSum = VectorXd::Zero(iterations);
  PeakSum2 = VectorXd::Zero(iterations);
	r=VectorXd::Zero(I);
	t=VectorXd::Zero(J);
	ax=VectorXd::Zero(J);
    VectorXd one=VectorXd::Ones(I);
    VectorXd zero1=VectorXd::Zero(I);
    VectorXd zero2=VectorXd::Zero(J);
    Z.conservativeResize(J,I);

	if (prompt == 'y' || prompt == 'Y')
	{
		for (int i = 0; i < Z.rows(); i++)
		{
			for (int j = 0; j < Z.cols(); j++)
			{
				Z(i,j) = s(i)*Z(i,j);
			}
		}
	}


	for (int i = 0; i < Z.rows(); i++)
	{
		norm_factor = 0;
		for (int j = 0; j < Z.cols(); j++)
		{
			if (j < threshold) {Z(i,j) = 0;}

			norm_factor += Z(i,j);
		}

		for (int j = 0; j < Z.cols(); j++)
		{
			if (norm_factor != 0) {
			Z(i,j) = Z(i,j)/norm_factor;
			}
			else Z(i,j) = 0;
		}
	}

    diff=VectorXd::Zero(I);
    R=Z.transpose();
    y.conservativeResize(I,1);
    sum=Update(Z*one,zero2,J);

	cout << "Starting MLEM..." << endl;

    /** ----------------------------------------------------
    * Maximum-Likelihood Expectation Maximization (MLEM)
    *   ----------------------------------------------------
    */

    for(k=0;k<iterations;k++){

	  x=x.cwiseProduct((Z*(y.cwiseQuotient(Update(R*x,zero1,I)))).cwiseQuotient(sum));

	  for(int i = threshold; i < J; i++) { t(i) = pow(x(i) - (x.sum()/(J-threshold)),2); }
	  dx =  dx + t;

	  r = (y - R*x);
	  residual(k) = sqrt(r.transpose()*r);

	  if (k%10==0) {cout << "Iteration(s): " << k << "\r" << flush;}

    if (iterations > 0 && y.sum() > 0)
      {
        diff = (R*x - diff);
        max_diff = diff.maxCoeff()/y.sum();
        //cout << k << " " << max_diff << endl;
        //if (max_diff < criteria) {break;}
        diff = R*x;

        // Max difference
        IndicatorMaxDiff(k) = max_diff;

        // Root-mean-square-error
        iter(k) = 0;
        for (int iii = 0; iii < x.size(); iii++)
        {
        iter(k) += (ToF(iii) - x(iii))*(ToF(iii) - x(iii));
        }
        iter(k) = sqrt(iter(k)/x.size());

        // Indicator function
        IndicatorFunc(k) = 0;
        float q, q_tot, tot_points;
        q_tot = 0; tot_points = 0;
        for (int iii = 0; iii < y.size(); iii++)
        {
          q = 0;
          for (int iiii = 0; iiii < x.size(); iiii++)
          {
            q += R(iii,iiii)*x(iiii);
          }
        IndicatorFunc(k) += (y(iii) - q)*(y(iii) - q);

        // Calc chi2
        if (y(iii) > 0) {
          chi2(k) += ((y(iii) - q)*(y(iii) - q))/(sqrt(y(iii))*sqrt(y(iii)));
        }
        tot_points += y(iii);
        q_tot += q;
        }
        chi2(k) /= tot_points;
        IndicatorFunc(k) = IndicatorFunc(k)/q_tot;


        PeakSum(k) = 0;
        PeakSum2(k) = 0;
        for (int n_peak = 56; n_peak <= 80; n_peak++)
        {
            PeakSum(k) += x(n_peak);
        }
        for (int n_peak = 28; n_peak <= 36; n_peak++)
        {
            PeakSum2(k) += x(n_peak);
        }




        //iter(k) = max_diff;
      }
    }

	dx = dx/(k+1);

    /** ----------------------------------------------------
    * Kevin's covariance matrix implementation
    *   ----------------------------------------------------
    */
    int I = min(sLength, lBins);
    int J = eBins;
    Z.conservativeResize(J,I);
    MatrixXd C = MatrixXd::Zero(J,J);   //Covariance matrix
    MatrixXd M_T = MatrixXd::Zero(I,J); //Transpose of Mji
    MatrixXd M = MatrixXd::Zero(J,I);   // Mji Matrix
    MatrixXd C_1 = MatrixXd::Zero(J,I); // Temporary Matrix
    //DumpMatrixToFile("./Z_2.dat", Z);

    for (int i = 0; i < Z.rows(); i++)
    {
        double norm_factor = 0;
        for (int j = 0; j < Z.cols(); j++)
        {
            if (j < threshold) {
                Z(i,j) = 0;
            }
            norm_factor += Z(i,j);
        }

        for (int j = 0; j < Z.cols(); j++)
        {
            if (norm_factor != 0)
            {
                Z(i,j) = Z(i,j) / norm_factor;
            }
            else
            {
                Z(i,j) = 0;
            }
        }
    }

    //DumpMatrixToFile("./Z_3.dat", Z);

    cout << "(Z) Inv. Response matrix has " << Z.cols() << " columns." << endl;
    cout << "(Z) Inv. Response matrix has " << Z.rows() << " rows." << endl;

    R = Z.transpose();


    //DumpMatrixToFile("./R_1.dat", R);

    cout << "(R) Response matrix has " << R.cols() << " columns." << endl;
    cout << "(R) Response matrix has " << R.rows() << " rows." << endl;


    //y.conservativeResize(I);

    cout << "y after resizing has " << y.cols() << " columns." << endl;
    cout << "y after resizing has " << y.rows() << " row." << endl;

    MatrixXd sum = Update(Z*one, zero2, J);
    VectorXd x_1 = VectorXd::Zero(J);
    MatrixXd x_11 = MatrixXd::Zero(J,J);
    MatrixXd y_1 = MatrixXd::Zero(J,J);
    VectorXd sum1 = Update(R*x, zero1, I);
    MatrixXd R_11 = MatrixXd::Zero(I,J);
    VectorXd std = VectorXd::Zero(J);
    MatrixXd std_op = MatrixXd::Zero(J,J);
    MatrixXd Corr = MatrixXd::Zero(J,J);  //Correlation Matrix

    /*-----------------------------------
    Covariance Matrix Calculation
    -------------------------------------------
    */

    x_1 = x.cwiseQuotient(sum); //x_1 vector stores the cwiseQuotient of x/(sum) here sum denotes efficiency
    x_11 = x_1.asDiagonal();   //Now I am converting the x_1 vector into a diagonal matrix of size J X J
    R_11 = R.cwiseQuotient(sum1.rowwise().replicate(J)); //In this step since sum1 = R*x is a vector of size I , I am replicating it J times in order to get the matrix of size I X J , so that I can take the cwiseQuotient with R of size (I X J)


    M_T = R_11*x_11; //Then multiplying by R_11*x_11 we will get the same answer as cwiseProduct

    y_1 = y.asDiagonal(); // // Converting y into a diagonal matrix of size J X J

    M = M_T.transpose();

    C_1 = M*y_1;

    C = C_1*M_T;  //Final Covariance Matrix


    /*----------------------------------
     Correlation Matrix
    --------------------------------------
     */
     std = C.diagonal().cwiseSqrt();
     std_op = std*std.transpose(); //Outer Product of standard deviation vector with it's transpose

     Corr = C.cwiseQuotient(std_op);

    //M = (y.transpose()).cwiseProduct(x);

    //((x.cwiseProduct(R)).transpose()).cwiseQuotient(R*x);

    //cout << M.rows() << " " << M.cols() << endl;

    //M.conservativeResize(M.rows(),M.rows());

    //V = y.asDiagonal();

    //C = M*V*M.transpose();

    out.open("covariance.out");
    out<<C;
    out.close();


    /** ----------------------------------------------------
    * Output the Unfolded Spectrum and errors
    *   ----------------------------------------------------
    */


    out.open("Unfolded.out");
    out<<x;
    out.close();

	  out.open("error.out");
    out<<dx;
    out.close();

	/** ----------------------------------------------------
    * Output the Estimate Spectrum
    *   ----------------------------------------------------
    */

	out.open("estimate.out");
	out<<R*x;
	out.close();

	/** ----------------------------------------------------
    * Output the Initial Spectrum
    *   ----------------------------------------------------
    */

	out.open("spectrum.in");
	out<<y;
	out.close();

  /** ----------------------------------------------------
    * Output the Peak Integral
    *   ----------------------------------------------------
    */

	out.open("peak_counts.out");
	out<<PeakSum;
	out.close();
  out.open("peak_counts2.out");
  out<<PeakSum2;
  out.close();

	/** ----------------------------------------------------
    * Output the residuals
    *   ----------------------------------------------------
    */

	out.open("residual.out");
	out<<residual;
	out.close();

  /** ----------------------------------------------------
    * Output the iterations and chi2
    *   ----------------------------------------------------
    */

  out.open("iterations_rmse.out");
  out<<iter;
  out.close();

  out.open("iterations_if.out");
  out<<IndicatorFunc;
  out.close();

  out.open("chi2.out");
  out<<chi2;
  out.close();

  out.open("maxdiff.out");
  out<<IndicatorMaxDiff;
  out.close();


    Stop = clock();
    cout << "Iterations run : " << k << endl;
    cout << "Elapsed time(s): " << (float)(Stop - Start)/CLOCKS_PER_SEC << endl;

	//cout << "Run again? (y/n)";
	//cin >> prompt;
  run = 0;
  //if(prompt == 'n') {run = 0;}

	}

  return 0;
}
