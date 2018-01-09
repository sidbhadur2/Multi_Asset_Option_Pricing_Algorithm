/*
 IE525 Numerical method in Finance
 Pricing Multi-asset option
 
 created by Joonha Yoon, Sid, Sunjie Hou
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>
#include <string.h>
#include <algorithm>
#include <iomanip>
#include "newmat.h"
#include "newmatio.h"
#include <chrono>
#include <random>
#include <time.h>
typedef std::chrono::high_resolution_clock Clock; // to use clock!
unsigned seed = (int)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
using namespace std;

vector <string> Ticker_name;
vector <double> Volatility;
vector <double> Initial_Price;
vector <double> Nu, up_factor, down_factor;
vector <double> temp;
double Risk_free_rate;
double Strike_price;
double Time_to_maturity;
double No_of_steps;
double **memo;
double temp_sum;

int No_of_Assets;
int No_of_trial;
Matrix Corr;
Matrix factor_index;
Matrix sign_assets;
Matrix Probability;
Matrix Unmapping;
double delta_T, discount;

int counter = 0;
int decoding_temp;
int mapping_value;
int unmapping_temp;
int prob_index;


double max(double a, double b) {
    return (b < a) ? a : b;
}

double even_odd_check(int k) { // change zeros to -1 for sign change
    if (k % 2 == 0) return 1;
    else return -1;
}

double summation(vector <double> k) { // vector sum
    temp_sum = 0;
    for (int i = 0; i < (int)k.size(); i++)
        temp_sum += k[i];
    return temp_sum;
}

double max_new_general(vector <double> gen_asset) { // find the highest stock price
    auto biggest = max_element(begin(gen_asset), end(gen_asset));
    return *biggest;
}




//////////////////////////////////// MEMOIZATION ////////////////////////////////////////
void make_memo() { // at kth step, the number of memo is k^n
    memo = new double *[(int)No_of_steps];
    for (int i = 0; i <= (int)No_of_steps; i++) { // make memo, first dimension refers a time step, and secon dimension refers a node number.
        memo[i] = new double[(int)pow((i + 1), (int)No_of_Assets)]; // as time step increases, the number of nodes increases in geometric process.
        for (int j = 0; j < pow((i + 1), (int)No_of_Assets); j++) {
            memo[i][j] = -1; // put -1 to all array which is invalid number for a option price.
        }
    }
}

int mapping(int time_step, Matrix index) { // Decoder for mapping N dimensional array to 1 D array.
    decoding_temp = 0;
    index = index + (time_step);
    for (int i = 1; i <= (int)No_of_Assets; i++) // here, we can get N decimal number. Each node has unique N decimal number.
        index(1, i) = index(1, i) / 2;
    
    for (int i = 0; i < (int)No_of_Assets; i++) // change this N decimal number to decimal number. this becomes the array index.
        decoding_temp += (int)index(1, (int)No_of_Assets - i)*(int)pow((time_step + 1), (int)No_of_Assets - i - 1);
    
    return decoding_temp;
}

double memoization_option_pricing(int k, Matrix index) { // k refers time step, and index refers the indice of upfactor of each stocks.
    
    if (memo[k][mapping(k, index)] != -1) 	return memo[k][mapping(k, index)]; // if there is a memoized value, return that value.
    
    if (k >= No_of_steps) { // if the algorithm reaches the end node, calculate the payoff.
        for (int i = 0; i < No_of_Assets; i++)
            temp[i] = Initial_Price[i] * pow(up_factor[i], index(1, i + 1));		// calculate the stock price at the end node.
        memo[k][mapping(k, index)] = max(max_new_general(temp) - Strike_price, 0.0); // Since we are pricing European max option, payoff is max(max(Si) - K, 0).
        return memo[k][mapping(k, index)];											// also, save the value at the corresponding memo.
    }
    else
    {
        memo[k][mapping(k, index)] = 0;
        for (int i = 0; i < pow(2, No_of_Assets); i++) // add up all path values. call the functions at the next time step (recursive algorithm)
            memo[k][mapping(k, index)] += discount*(Probability(i + 1, 1)*memoization_option_pricing(k + 1, (sign_assets.Row(i + 1) + index)));
        
        return memo[k][mapping(k, index)];
    }
}



///////////////////////////////DYNAMIC PROGRAMMING///////////////////////////////////////

void clear_memo() { // clear memo for reusing it in dynamic programming,
    for (int i = 0; i <= (int)No_of_steps; i++)
        for (int j = 0; j < pow((i + 1), (int)No_of_Assets); j++)  // put zeros to all.
            memo[i][j] = 0;
}

Matrix unmapping(int mapping_index, int time_step) { // this encoder change decimal numbers to the N decimal number.
    // and generates corresponding index matrix.
    for (int i = 0; i < No_of_Assets; i++) {
        unmapping_temp = (int)pow((time_step + 1), (int)No_of_Assets - i - 1);
        Unmapping(1, No_of_Assets - i) = mapping_index / unmapping_temp; // change decimal number to N decimal number.
        mapping_index = mapping_index%unmapping_temp; //  for example, decimal number 6 is 1*2^2 + 1*2^1 + 0*2^0 which is 110 in binary number.
    }
    Unmapping = Unmapping * 2 - time_step; // multiply 2 and subtract time_step to make it balanced index number.
    
    return Unmapping;
}

void calculating_end_val() {
    Matrix index(1, No_of_Assets);
    
    for (int i = 0; i <  pow((No_of_steps + 1), (int)No_of_Assets); i++) { // the number of end nodes is pow((No_of_steps + 1), (int)No_of_Assets)
        memo[(int)No_of_steps][i] = 0;
        index = unmapping(i, (int)No_of_steps);		// calculate the up_factors' indice
        
        for (int j = 0; j < No_of_Assets; j++)
            temp[j] = Initial_Price[j] * pow(up_factor[j], index(1, j + 1));
        
        memo[(int)No_of_steps][i] = max(max_new_general(temp) - Strike_price, 0.0); // calculate the payoff and save it
    }
}

double dynamic_option_pricing() {
    Matrix index(1, No_of_Assets);
    
    for (int i = (int)No_of_steps-1; i >= 0; i--) { // algorithm moves backward.
        for (int j = 0; j < pow((i + 1), (int)No_of_Assets); j++) { // search all nodes in the time step
            for (int k = 0; k < (int)pow(2, No_of_Assets); k++) { // add up all branches
                index = unmapping(j, i); // change the array index to the node's upfactor index.
                memo[i][j] += discount*(Probability(k + 1, 1))*memo[i+1][mapping(i + 1, index + sign_assets.Row(k + 1))];
            }	//find all branches of the node by considering all possible index change. and add up all (values *  risk_neutral probability)
        }
    }
    return memo[0][0];
}

////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////// MONTE CARLO SIMULATION WITH LATTICE MODEL ///////////////////////////

double get_uniform() // get Ui ~ U[0,1]
{
    std::uniform_real_distribution <double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

int pick_path(double prob) { //  prob ~ U[0,1]
    prob_index = 1;
    for (int i = 1; i <= (int)pow(2, No_of_Assets); i++) {
        if ((prob - Probability(i, 1) < 0)) break; // divide [0,1] into N with N many risk neutral probabilities.
        else prob = prob - Probability(i, 1);	//find corresponding probability.
        prob_index++;
    }
    return prob_index;
}

double Monte_carlo_option_pricing() {
    
    Matrix index(1, No_of_Assets);
    int prob_index;
    double monte_result = 0;
    
    for (int trial = 0; trial < No_of_trial; trial++) {
        index = factor_index;
        for (int step = 0; step < No_of_steps; step++) {
            prob_index = pick_path(get_uniform()); // randomly pick a path. the probability of picking up a path is the same with the corresponding risk neutral probability.
            index = sign_assets.Row(prob_index) + index; // change the up factor index
        }
        
        for (int i = 0; i < No_of_Assets; i++)
            temp[i] = Initial_Price[i] * pow(up_factor[i], index(1, i + 1)); // calculate the payoff
        
        monte_result += exp(-Risk_free_rate * Time_to_maturity)* max(max_new_general(temp) - Strike_price, 0.0);
        
    }
    
    return (monte_result/(double)No_of_trial);
}


///////////////////// MONTE CARLO SIMULATION WITH LATTICE MODEL ///////////////////////////

void Get_Parameter_data() {
    ifstream file_1("/Users/siddharthbhaduri/Desktop/Work/Spring-2017/IE-525-Numerical_Methods/Multi_Asset_Option_Project/info.csv");
    string value;
    bool read = true;
    bool pass_column_names = true;
    while (file_1.good()) // read csv file
    {
        if (pass_column_names) {
            for (int j = 0; j < 6; j++)
                getline(file_1, value, ',');
            pass_column_names = false;
        }
        getline(file_1, value, ',');
        getline(file_1, value, ',');
        if (string(value) == "") break; // if Ticker name is empty, stop reading
        Ticker_name.push_back(string(value));//read Ticker name
        getline(file_1, value, ',');
        Initial_Price.push_back(atof(string(value).c_str())); // read initial prices
        getline(file_1, value, ',');
        Volatility.push_back(atof(string(value).c_str())); // read volatility
        getline(file_1, value, ',');
        if (read) 	Risk_free_rate = atof(string(value).c_str()); // read risk_free rate
        getline(file_1, value, ',');
        if (read) 	Time_to_maturity = atof(string(value).c_str());
        getline(file_1, value);
        if (read) {
            Strike_price = atof(string(value).c_str());
            read = false;
        }
    }
    
    No_of_Assets = (int)Ticker_name.size(); // calculate the number of assets N
    No_of_Assets = 4;
    Matrix A(No_of_Assets, No_of_Assets); // make NxN array
    
    ifstream file_2("/Users/siddharthbhaduri/Desktop/Work/Spring-2017/IE-525-Numerical_Methods/Multi_Asset_Option_Project/corr.csv");
    if (file_2.good())
    {
        for (int i = 1; i <= No_of_Assets + 1; i++) // get correlation coefficient data.
            for (int j = 1; j <= No_of_Assets+1; j++) {
                if (j == No_of_Assets + 1) getline(file_2, value);
                else getline(file_2, value, ',');
                if ((i >= 2)&&(j>=2)) {
                    A(i - 1, j - 1) = atof(string(value).c_str());
                }
            }
        Corr = A;
    }
    cout << "#####################" << endl;
    cout << " Data Read Complete " << endl;
    cout << "#####################" << endl;
}

void Print_data() {
    cout << endl << "Ticker	S_0	Rf	Sigma	T	K" << endl;
    for (int i = 0; i < (int)Ticker_name.size(); i++) {
        cout << Ticker_name[i] << "	" << Initial_Price[i] << "	" << Risk_free_rate << "	" << Volatility[i] << "	" << Time_to_maturity << "	" << Strike_price << endl;
    }
    cout << endl << "Correlation Coefficient Matrix :" << endl;
    cout << endl << setw(10) << setprecision(5) << Corr;
}


void Initializing() {
    
    int No_of_prob = (int)pow(2, No_of_Assets); // the number of risk neutral probabilities.
    int No_of_corr = No_of_Assets*(No_of_Assets - 1) / 2; // the number of rho
    
    Matrix sign_rho(No_of_prob, No_of_corr);
    Matrix Nu_devided_by_sigma(No_of_Assets, 1);
    Matrix Corr_column_vector(No_of_corr, 1);
    Matrix A(No_of_prob, No_of_Assets);
    Matrix B(No_of_prob, 1);
    Matrix C(1, No_of_Assets);
    sign_assets = A;
    Probability = B;
    factor_index = C;
    Unmapping = C;
    
    delta_T = Time_to_maturity / No_of_steps;
    discount = exp(-Risk_free_rate*delta_T);
    
    for (int i = 0; i < No_of_Assets; i++) {
        Nu.push_back(Risk_free_rate - (0.5)*pow(Volatility[i], 2));
        up_factor.push_back(exp(Volatility[i] * sqrt(delta_T)));
        down_factor.push_back(pow(up_factor[i], -1));
        factor_index(1, i + 1) = 0;
        temp.push_back(0);
    }
    
    for (int i = 1; i <= No_of_Assets; i++)
        Nu_devided_by_sigma(i, 1) = Nu[i - 1] / Volatility[i - 1]; // calculate (r-sigma^2/2)/sigma
    
    int count = 1;
    for (int j = 1; j <= No_of_Assets - 1; j++)
        for (int z = j + 1; z <= No_of_Assets; z++) { // line up all correlation coefficients for the matrix calculation
            Corr_column_vector(count, 1) = Corr(j, z);
            count++;
        }
    
    int temp;
    for (int i = 1; i <= No_of_Assets; i++)
        for (int j = 1; j <= No_of_prob; j++) { // make a matrix with all possible sign changes for Nu/sigma
            temp = (j - 1) / ((int)pow(2, No_of_Assets - i));
            sign_assets(j, i) = even_odd_check(temp);
        }
    
    for (int i = 1; i <= No_of_prob; i++) {
        count = 1;
        for (int j = 1; j <= No_of_Assets - 1; j++) // make a matrix with all corresponding sign changes for rho
            for (int z = j + 1; z <= No_of_Assets; z++) {
                sign_rho(i, count) = sign_assets(i, j) * sign_assets(i, z);
                count++;
            }
    }
    
    Probability = (1 / (pow(2, No_of_Assets)) * (1 + (sign_assets*Nu_devided_by_sigma) *sqrt(delta_T) + sign_rho * Corr_column_vector));
}

///////////////////// RECURSIVE ALGORITM ///////////////////////////

double binomial_option_pricing(int k, Matrix index){	// memory problem !
    if (k >= No_of_steps) {
        for (int i = 0; i < No_of_Assets; i++)
            temp[i] = Initial_Price[i] * pow(up_factor[i], index(1,i+1));
        return max(max_new_general(temp) - Strike_price, 0.0);
    }
    else
    {
        double total = 0;
        for (int i = 0; i < pow(2, No_of_Assets); i++)
            total += discount*(Probability(i + 1, 1)*binomial_option_pricing(k + 1, (sign_assets.Row(i + 1) + index)));
        return total;
    }
}
/////////////////////////////////////////////////////////////////////////////


int main() {
    
    ////////for testing////////////////
    Risk_free_rate = 0.0213;
    No_of_steps = 7;
    No_of_trial = 10000;
    //Strike_price = 100;
    
    /////////////////////////////////////
    
    cout << "---------------------------------------" << endl;
    Get_Parameter_data();
    Initializing();
    Print_data();
    cout << "---------------------------------------" << endl << endl;
    
    
    make_memo();
    auto t_1 = Clock::now();
    cout << "Option price using Recursion : " << binomial_option_pricing(0,factor_index);
    auto t_2 = Clock::now();
    cout << "		" << chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9) << " sec"<< endl;
    cout << endl;
   
    /*
    t_1 = Clock::now();
    cout << "Option price using Memoization :         " << memoization_option_pricing(0, factor_index);
    t_2 = Clock::now();
    cout << "		" << chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9) << " sec" << endl;
    cout << endl;
    
    
    t_1 = Clock::now();
    cout << "Option price using Monte Carlo :         "<< Monte_carlo_option_pricing();
    t_2 = Clock::now();
    cout << "		" << chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9) << " sec" << endl;
    cout << endl;
    
    clear_memo();
    calculating_end_val();
    t_1 = Clock::now();
    cout << "Option price using Dynamic Programming : "<< dynamic_option_pricing();
    t_2 = Clock::now();
    cout << "		" << chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9) << " sec" << endl;
    
    */
    return 0;
    
}


