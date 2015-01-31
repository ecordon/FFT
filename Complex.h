// Complex.h
// Complex numbers class
// Edwin Cordon

#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>

class Complex{

	public: 
		
        // Constructors  

        Complex(double magnitude, double phase);
        Complex();
			
		// Overloaded operators with Complex types

        Complex operator+(const Complex &) const;
        Complex operator-(const Complex &) const;
        Complex operator*(const Complex &) const;
        Complex operator/(const Complex &) const;
        
        Complex &operator+=(const Complex &);
		Complex &operator-=(const Complex &);
		Complex &operator*=(const Complex &);
		Complex &operator/=(const Complex &);

        Complex &operator=(const Complex &);

		// Overloaded operators with double types
        
        Complex operator+(const double z) const;
        Complex operator*(const double z) const;
		Complex &operator=(const double z);


        // Methods

		void setPolar(double magnitude, double phase);
		double getMagnitude();
		double getPhase();

		void multiplyBy(double a);
		void shiftBy(double theta);

		void setRect(double real, double imaginary);
		double getReal();
		double getImaginary();

        void printPolar();
        void printRectangular();

	private:
		
        double Re;
		double Im;
		double magnitude;
		double phase;
};

#endif

