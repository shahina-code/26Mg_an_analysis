/*************************************************************************
 *
 *  Filename: UncertainNumber.h
 *
 *  Description:
 *
 *
 *	  - This was adpated from a java class from Dale Visser (Yale).
 *
 *	Author(s):
 *     Michael T. Febbraro
 *
 *  Creation Date: 7/2/2016
 *  Last modified: 3/31/2019
 *
 * -----------------------------------------------------
 * 	Physics Divsion - Oak Ridge National Laboratory
 *	Oak Ridge, TN, USA
 *
 */
#include <iostream>
#include <fstream>
#include <ctime>
#include <signal.h>
#include <string>
#include <cstring>
#include <math.h>

using namespace std;

class UncertainNumber {
private :

protected :

public :

	double	value;
	double	error;

	UncertainNumber(double val = 0, double er = 0) : value(val) , error(er) {};

	~UncertainNumber() {};

	// Standard functions
	UncertainNumber operator+ (const UncertainNumber & x) const
		{return UncertainNumber(value + x.value, sqrt(error * error + x.error * x.error));}
	UncertainNumber operator+ (const double & x) const
		{return UncertainNumber(value + x, error);}
	UncertainNumber operator- (const UncertainNumber & x) const
		{return UncertainNumber(value - x.value, sqrt(error * error + x.error * x.error));}
	UncertainNumber operator- (const double & x) const
		{return UncertainNumber(value - x, error);}
	UncertainNumber operator* (const UncertainNumber & x) const
		{return UncertainNumber(value * x.value,
		sqrt(error * error * x.value * x.value + x.error * x.error * value * value));}
	UncertainNumber operator* (const double x) const
		{return UncertainNumber(value*x, sqrt(x*error));}
	UncertainNumber operator/ (const UncertainNumber & x) const
		{return UncertainNumber(value / x.value, sqrt(error * error * x.value * x.value + x.error * x.error * value * value));}

	// Power functions
	//UncertainNumber operator^ (const UncertainNumber & x) const
	//	{return UncertainNumber(pow(value, x.value), pow(value, x.value)*sqrt(pow(((x.value/value)*error), 2) + pow((log(value)*x.error)), 2);}
	UncertainNumber operator^ (const double & x) const
		{return UncertainNumber(pow(value, x), ((pow(value, x)*x*error)/value));}

	// Need Trig functions

	UncertainNumber ln() const {return UncertainNumber(log(value), error/value);}

	double getValue() const {return value;}
	double getError() const {return error;}

	void set(double val, double er) {value = val; error = er;}

	string toString() { return to_string(value) + "(" + to_string(error) + ")"; }
};
