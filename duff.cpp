#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using force_function = std::vector<double> (*)(const std::vector<double>&, const std::vector<double>&, double, double, double, double, double, double); // e' pensata per la forza dell'oscillatore di duffing

void print_all_i_need(const std::vector<double>& x, const std::vector<double>& v, const double t) { // per stampare posizione, velocita' e tempo
    std::cout << "x: ";
    for (double element : x) {
        std::cout << element << " ";
    }
    std::cout << "v: ";
    for (double element : v) {
        std::cout << element << " ";
    }
    std::cout << "t: " << t << "\n";
}

std::vector<double> prod(const std::vector<double>& x, const double s) { // moltiplicazione di un vettore per uno scalare

    std::vector<double> product;

    for (double element : x) {
        product.push_back(element * s);
    }

    return product;
}

std::vector<double> sum(const std::vector<double>& x, const std::vector<double>& y) { // somma tra due vettori

    std::vector<double> sum;

    // controllo se x e y hanno la stessa dimensione
    if (x.size() != y.size()) {
        // errore se hanno dimensioni diverse
        std::cerr << "Error: Vectors x and y must have the same size.\n";
        // Restituisci un vettore vuoto
        return sum;
    }

    for (size_t i = 0; i < x.size(); i++) {
        sum.push_back(x[i]+y[i]);
    }

    return sum;
}

std::vector<double> duff_force(const std::vector<double>& x, const std::vector<double>& v, const double t, const double a, const double b, const double g, const double d, const double o) { // dall'equazione di duffing

    std::vector<double> duffing;

    for (size_t i = 0; i < x.size(); i++) {
        double duffing_i = d*v[i] + a*x[i] + b*x[i]*x[i]*x[i] - g*std::cos(o*t);
        duffing.push_back(duffing_i);
    }

    return duffing;
}

void euler(std::vector<double>& x, std::vector<double>& v, force_function force, double dt, double t_max, const double a, const double b, const double g, const double d, const double o) {

    int n_step = t_max / dt;
    double t = 0;
    for (int i=0; i < n_step; i++) {
        t += dt;
        x = sum(x, prod(v,dt));
        v = sum(v, prod(duff_force(x, v, t, a, b, g, d, o), -dt));
        if (n_step % 10000 == 0){print_all_i_need(x,v,t);} // stampo ogni 10000 passi temporali
    }
}

int main (int argc, char* argv[]){ //prendo i dati iniziali da terminale inserendoli dopo il lancio del programma nella stessa riga e.g.: ./euler 1.0 0.0 0.0 0.0 1.4142 0.0

   if (argc < 3 || (argc % 2) != 1) { // controllo di avere un numero pari di parametri
        std::cerr << "Usage: " << argv[0] << " value1_1 value1_2 ... valueN-1_1 valueN-1_2\n";
        return 1;  // ritorna con errore
    }

    // ottengo il numero di elementi
    int num_elements = (argc - 1) / 2;

    // dichiaro i vettori x e v
    std::vector<double> x_0(num_elements);
    std::vector<double> v_0(num_elements);

    // popolo i vettori x e v dagli argomenti dopo il nome del programma
    for (int i = 0; i < num_elements; ++i) {
        std::istringstream(argv[1 + i]) >> x_0[i];
        std::istringstream(argv[1 + num_elements + i]) >> v_0[i];
    }

    //std::vector<double> x_0 = {1.0, 0.0, 0.0};
    //std::vector<double> v_0 = {0.0, 1.4142, 0.0};
    //std::vector<double> x;
    //std::vector<double> v;
    // parametri del modello e dati iniziali
    double a = -1;
    double b = 0.25;
    double g = 2.5;
    double d = 0.1;
    double o = 2;
    double dt = 0.00001;
    double t = 0;

    print_all_i_need(x_0, v_0, t); // stampo i dati iniziali
    euler(x_0, v_0, duff_force, dt, 100, a, b, g, d, o); // seguo il punto materiale per 100 unita' di tempo
}
