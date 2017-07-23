#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "Eigen\Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// TODO - complete this function
vector<double> JMT(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficient in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */
    
    MatrixXd A = MatrixXd(3, 3);
	VectorXd B = MatrixXd(3,1);	

	double t2 = t * t;
	double t3 = t * t2;
	double t4 = t * t3;
    double t5 = t * t4;

    A <<   t3,    t4,    t5,
         3*t2,  4*t3,  5*t4,
         6*t,  12*t2, 20*t3;
		
	    
    B << end[0] - (start[0] + start[1] * t + 0.5 * start[2] * t2),
         end[1] - (start[1] + start[2] * t),
         end[2] - start[2];
			    
	
	VectorXd C = A.inverse() * B;
	
	double a0 = start[0];
    double a1 = start[1];
    double a2 = start[2] * 0.5;
    double a3 = C[0];
    double a4 = C[1];
    double a5 = C[2];

	vector<double> result = {a0, a1, a2, a3, a4, a5}; 
	
    return result;
    
}

bool close_enough(vector< double > poly, vector<double> target_poly, double eps=0.01) {


	if(poly.size() != target_poly.size())
	{
		cout << "your solution didn't have the correct number of terms" << endl;
		return false;
	}
	for(int i = 0; i < poly.size(); i++)
	{
		double diff = poly[i]-target_poly[i];
		if(abs(diff) > eps)
		{
			cout << "at least one of your terms differed from target by more than " << eps << endl;
			return false;
		}

	}
	return true;
}
	
struct test_case {
	
		vector<double> start;
		vector<double> end;
		double T;
};

vector< vector<double> > answers = {{0.0, 10.0, 0.0, 0.0, 0.0, 0.0},{0.0,10.0,0.0,0.0,-0.625,0.3125},{5.0,10.0,1.0,-3.0,0.64,-0.0432}};

int main() {

	//create test cases

	vector< test_case > tc;

	test_case tc1;
	tc1.start = {0,10,0};
	tc1.end = {10,10,0};
	tc1.T = 1;
	tc.push_back(tc1);

	test_case tc2;
	tc2.start = {0,10,0};
	tc2.end = {20,15,20};
	tc2.T = 2;
	tc.push_back(tc2);

	test_case tc3;
	tc3.start = {5,10,2};
	tc3.end = {-30,-20,-4};
	tc3.T = 5;
	tc.push_back(tc3);

	bool total_correct = true;
	for(int i = 0; i < tc.size(); i++)
	{
		vector< double > jmt = JMT(tc[i].start, tc[i].end, tc[i].T);
		bool correct = close_enough(jmt,answers[i]);
		total_correct &= correct;

	}
	if(!total_correct)
	{
		cout << "Try again!" << endl;
	}
	else
	{
		cout << "Nice work!" << endl;
	}

	return 0;
}