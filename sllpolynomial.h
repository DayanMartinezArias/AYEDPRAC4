// AUTOR: 
// FECHA: 
// EMAIL: 
// VERSION: 2.0
// ASIGNATURA: Algoritmos y Estructuras de Datos
// PRÁCTICA Nº: 4
// ESTILO: Google C++ Style Guide
// COMENTARIOS:
// 

#ifndef SLLPOLYNOMIAL_H_
#define SLLPOLYNOMIAL_H_

#include <iostream>
#include <math.h>  // fabs, pow


#include "pair_t.h"
#include "sll_t.h"
#include "vector_t.h"


#define EPS 1.0e-6

typedef pair_t<double> pair_double_t;  // Campo data_ de SllPolynomial
typedef sll_node_t<pair_double_t> SllPolyNode;  // Nodos de SllPolynomial

// Clase para polinomios basados en listas simples de pares
class SllPolynomial : public sll_t<pair_double_t> {
 public:
  // constructores
  SllPolynomial(void) : sll_t() {};
  SllPolynomial(const vector_t<double>&, const double = EPS);

  // destructor
  ~SllPolynomial() {};

  // E/S
  void Write(std::ostream& = std::cout) const;
  
  // operaciones
  double Eval(const double) const;
  bool IsEqual(const SllPolynomial&, const double = EPS) const;
  void Sum(const SllPolynomial&, SllPolynomial&, const double = EPS);
  void lastElem();
  void sumFirstTwo();
  void writeInx(double val);
  void mean();
  void max();
  void count_mon();
  void consec();
  void sumaPares();
  void printZero();
  void Penultimo();
  void val_medio();
  void multiplosTres();
  void Integral(double lim_inf, double lim_sup);
  void Inversa();
};


bool IsNotZero(const double val, const double eps = EPS) {
  return fabs(val) > eps;
}

// FASE II
// constructor
SllPolynomial::SllPolynomial(const vector_t<double>& v, const double eps) {
for(int i = v.get_size() - 1; i >= 0; i--) {
    if(IsNotZero(v.at(i))) {
      pair_double_t monomio(v.at(i), i);
      SllPolyNode* nodo = new SllPolyNode(monomio);
      this->push_front(nodo);
    }
  }

}

// E/S
void SllPolynomial::Write(std::ostream& os) const {
  os << "[ ";
  bool first{true};
  SllPolyNode* aux{get_head()};
  while (aux != NULL) {
    int inx{aux->get_data().get_inx()};
    double val{aux->get_data().get_val()};
    if (val > 0)
      os << (!first ? " + " : "") << val;
    else
      os << (!first ? " - " : "-") << fabs(val);
    os << (inx > 1 ? " x^" : (inx == 1) ? " x" : "");
    if (inx > 1)
      os << inx;
    first = false;
    aux = aux->get_next();
  }
  os << " ]" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SllPolynomial& p) {
  p.Write(os);
  return os;
}


// Operaciones con polinomios

// FASE III
// Evaluación de un polinomio representado por lista simple
double SllPolynomial::Eval(const double x) const {
  double result{0.0};
  SllPolyNode* temp = this->get_head();
  while(temp != NULL) {
    result += temp->get_data().get_val() * pow(x, temp->get_data().get_inx());
    temp = temp->get_next();
  }

  return result;
}

// Comparación si son iguales dos polinomios representados por listas simples
bool SllPolynomial::IsEqual(const SllPolynomial& sllpol,
			    const double eps) const {
  bool differents = false;
  SllPolyNode* temp = this->get_head();
  SllPolyNode* temp2 = sllpol.get_head();
  while(temp != NULL) {
    if(fabs(((temp->get_data().get_val()) - (temp2->get_data().get_val())) > eps) || ((temp->get_data().get_inx()) != (temp2->get_data().get_inx()))){
      return differents;
    }
    temp = temp->get_next();
    temp2 = temp2->get_next();
    if((temp2 == NULL && temp != NULL) || (temp2 != NULL && temp == NULL)) {
      return differents;
    }
  }
  return !differents;
}

// FASE IV
// Generar nuevo polinomio suma del polinomio invocante mas otro polinomio
void SllPolynomial::Sum(const SllPolynomial& sllpol, 
                       SllPolynomial& sllpolsum, 
                       const double eps) {
 
  SllPolynomial auxSllPolSum;                            
  SllPolyNode* aux = get_head();                         
  SllPolyNode* aux2 = sllpol.get_head();               
  pair_double_t pair;                                     
  double sum = 0.0;

  // si tienen los mismos exponentes se suman
  while (aux != NULL || aux2 != NULL) {                   
    if (aux != NULL && aux2 != NULL) {                   
      if (aux->get_data().get_inx() == aux2->get_data().get_inx()) { 
        sum = aux->get_data().get_val() + aux2->get_data().get_val();
        if (IsNotZero(sum, eps)) {                        
          pair.set(sum, aux->get_data().get_inx());    
          auxSllPolSum.push_front(new SllPolyNode(pair));
        }
        if (aux != NULL) aux = aux->get_next();         
        if (aux2 != NULL) aux2 = aux2->get_next();     
      } else if (aux->get_data().get_inx() > aux2->get_data().get_inx()) { 
        sum = aux2->get_data().get_val();             
        if (IsNotZero(sum, eps)) {                      
          pair.set(sum, aux2->get_data().get_inx()); 
          auxSllPolSum.push_front(new SllPolyNode(pair));
        }
        if (aux2 != NULL) aux2 = aux2->get_next();
      } else {
        sum = aux->get_data().get_val();              
        if (IsNotZero(sum, eps)) {                       
          pair.set(sum, aux->get_data().get_inx());    
          auxSllPolSum.push_front(new SllPolyNode(pair));
        }

        if (aux != NULL) aux = aux->get_next();         
      }
                                                     
    } else if (aux == NULL) {                           
      sum = aux2->get_data().get_val();                
      if (IsNotZero(sum, eps)) {                        
        pair.set(sum, aux2->get_data().get_inx());    
        auxSllPolSum.push_front(new SllPolyNode(pair)); 
      }
      if (aux2 != NULL) aux2 = aux2->get_next();         

    } else if (aux2 == NULL) {                          
      sum = aux->get_data().get_val();                   
      if (IsNotZero(sum, eps)) {                        
        pair.set(sum, aux2->get_data().get_inx());      
        auxSllPolSum.push_front(new SllPolyNode(pair)); 
      }

      if (aux != NULL) aux = aux->get_next();           
    }
  }

  while (!auxSllPolSum.empty()) {                      
    sllpolsum.push_front(auxSllPolSum.pop_front());     
  }
}

//Hecho en casa para practicar

//Mostrar último elemento
void SllPolynomial::lastElem() {
  SllPolyNode* temp = this->get_head();
  while(temp->get_next() != NULL) {
    temp = temp->get_next();
  }
  std::cout << temp->get_data().get_val() << std::endl;
}

//Sumar los dos primeros
void SllPolynomial::sumFirstTwo() {
  double sum = 0.0;
  SllPolyNode* temp = this->get_head();
  for(int i = 0; i < 2; ++i) {
    sum += temp->get_data().get_val();
    temp = temp->get_next();
  }
  std::cout << sum << std::endl;
}


//Indice
void SllPolynomial::writeInx(double val) {
  SllPolyNode* temp = this->get_head();
  while(temp->get_next() != NULL) {
    if(val == temp->get_data().get_val()) {
      std::cout << temp->get_data().get_inx() << std::endl;
      return;
    }
    temp = temp->get_next();
  }
  std::cout << "No encontrado" << std::endl;
}

//media
void SllPolynomial::mean() {
  double sum = 0.0;
  double counter = 0.0;
  SllPolyNode* temp = this->get_head();
  while(temp != NULL) {
    sum += temp->get_data().get_val();
    temp = temp->get_next();
    counter += 1.0;
  }
  double mean = sum / counter;
  std::cout << mean << std::endl;
}

//contar monomio
void SllPolynomial::count_mon() {
  double counter = 0.0;
  SllPolyNode* temp = this->get_head();
  while(temp != NULL) {
    temp = temp->get_next();
    counter += 1.0;
  }
  std::cout << counter << std::endl;
}

//Coeficiente maximo
void SllPolynomial::max() {
  double max = 0.0;
  SllPolyNode* temp = this->get_head();
  while(temp != NULL) {
    if(temp->get_data().get_val() > max) {
      max = temp->get_data().get_val();
    }
    temp = temp->get_next();
  }
}

//suma pares
void SllPolynomial::sumaPares() {
  double sum = 0.0;
  SllPolyNode* temp = this->get_head();
  while(temp != NULL) {
    if(temp->get_data().get_inx() % 2 == 0) {
      sum += temp->get_data().get_val();
    }
    temp = temp->get_next();
  }
}

// Imprimir los consecutivos
void SllPolynomial::consec() {
  SllPolyNode* temp = this->get_head();
  while(temp->get_next() != NULL) {
    if(fabs(temp->get_next()->get_data().get_inx() - temp->get_data().get_inx()) == 1) {
      std::cout << "(" << temp->get_data().get_val() << ", " << temp->get_next()->get_data().get_val() << ") " << std::endl;
    }
    temp = temp->get_next();
  }
}

//Imprimir los 0s
void SllPolynomial::printZero() {
  SllPolyNode* temp = this->get_head();
  while(temp->get_next() != NULL) {
    int i = temp->get_next()->get_data().get_inx() - temp->get_data().get_inx();
    if(i > 1) {
      for(int j = 1; j < i; ++j) {
        std::cout << "0x^" << (temp->get_data().get_inx() + j) << " ";
      }
    }
    temp = temp->get_next();
  }
  std::cout << std::endl;
}

//Imprimir penultimo elemento
void SllPolynomial::Penultimo() {
  SllPolyNode* temp = this->get_head();
  while(temp->get_next()->get_next() != NULL) {
    temp = temp->get_next();
  }
  std::cout << temp->get_data().get_val() << std::endl;
}

//Coeficiente que esta en el medio
void SllPolynomial::val_medio() {
  int counter = 0;
  SllPolyNode* temp = this->get_head();
  while(temp != NULL) {
    temp = temp->get_next();
    counter += 1;
  }
  temp = this->get_head();
  if(counter % 2 != 0) {
    for(int i = 0; i < counter/2; ++i) {
      temp = temp->get_next();
    }
  } else {
    std::cout << "El polinomio tiene un número par de monomios" << std::endl;
    return;
  }
  std::cout << temp->get_data().get_val() << std::endl;
}

//Sumar o imprimir todos los multiplos de x numero, en este caso sumar, pero imprimir es casi lo mismo cambiando unas cositas
void SllPolynomial::multiplosTres() {
  double sum = 0.0;
  SllPolyNode* temp= this->get_head();
  while(temp != NULL) {
    if(temp->get_data().get_inx() % 3 == 0) {
      sum += temp->get_data().get_val();
    }
    temp = temp->get_next();
  }
  std::cout << sum << std::endl;
}

//Supongo que ta bien, I guess
void SllPolynomial::Integral(double lim_inf, double lim_sup) {
  double sum = 0.0;
  SllPolyNode* aux = this->get_head();
  std::cout << "PRIMITIVA: " << std::endl;
  while(aux != NULL) {
    double new_inx = static_cast<double>(aux->get_data().get_inx() + 1);
    double new_coef = (aux->get_data().get_val()) / (new_inx);
    if(aux->get_next() == NULL) {
      std::cout << new_coef << "X^" << new_inx << std::endl;
    } else {
      std::cout << new_coef << "X^" << new_inx << " + ";
    }
    sum += (pow(lim_sup, new_inx) * new_coef) - (pow(lim_inf, new_inx) * new_coef);
    aux = aux->get_next();
  }
  std::cout << "INTEGRAL DEFNIDA :" << std::endl;
  std::cout << sum << std::endl;
}

//recorrido inversa


#endif  // SLLPOLYNOMIAL_H_
