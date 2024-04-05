#include<iostream>

double pot(double number, int exp);

int main() {
  std::cout << pot(2, 4) << std::endl;
}
 
double pot(double number, int exp) {
  if(exp == 1 ) {
    return number;
  }
  return number * pot(number, exp - 1);
}

