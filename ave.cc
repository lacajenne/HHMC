#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

void stat( const vector<double>& data, double& ave, double& error );
void jack( const vector<double>& data, double& ave, double& error );

int main(int argc, char *argv[]){

  int task;

  if(argc != 4){
    cerr << "Usage: ncols file task(0:stat, 1:jack) \n";
    return 1;
  }
	
  int ncols = atoi(argv[1]);
  ifstream data_file(argv[2]);
  if(!data_file || ncols <= 0){
    cerr << "Wrong parameters.\n";
    return 1;
  }
  task = atoi(argv[3]);
  
  vector<double>* data = new vector<double>[ncols];
  double tmp;
	
  int flag = 1;
  while(flag){
    for(int i = 0; i < ncols; i++){
      data_file >> tmp;
      if(data_file.eof()){ flag = 0; break; }
      data[i].push_back(tmp);
    }
  }
  
  double average;
  double error;
  for(int i = 0; i < ncols; i++){
    if(task==0){
      stat(data[i], average, error);
    } else if(task==1){
      jack(data[i], average, error);
    } else {
      cerr << "WRONG CHOICE" << endl;
    }
    cout << average << " " << error << "\n";
  }
  
  data_file.close();
  
  cerr << "ave on " << data[0].size() << "\n";
  
  return 0;
}


void stat( const vector<double>& data, double& ave, double& error )
{
  double* bin_ave;
  unsigned long int ndata = data.size();
  bin_ave = new double[ndata];
  
  ave = 0.0;
  for(unsigned long int n = 0; n < ndata; n++) ave += data[n];
  ave /= ndata;
	
  error = 0.0;
  for(unsigned long int bin_size = 1; bin_size <= ndata; bin_size ++){
    int nbins = ndata/bin_size;
    for(unsigned long int i = 0; i < nbins; i++){
      bin_ave[i] = 0.0;
      for(unsigned long int k = 0; k < bin_size; k++)
	bin_ave[i] += data[k+i*bin_size];
      bin_ave[i] /= bin_size;
    }
    double tmp = 0.0;
    for(unsigned long int i = 0; i < nbins; i++)
      tmp += ( bin_ave[i] - ave ) * ( bin_ave[i] - ave );
    tmp = sqrt(tmp)/nbins;
    if(tmp > error) error = tmp;
  }
  
  delete[] bin_ave;
}


void jack(const vector<double>& data, double& ave, double& error)
{
  double* bin_ave;
  unsigned long int ndata = data.size();
  bin_ave = new double[ndata];
  
  ave = 0.0;
  for(unsigned long int n = 0; n < ndata; n++) ave += data[n];
  ave /= ndata;  

  for(int n=0; n<ndata; n++){
    bin_ave[n] = 0.0;
    for(int k=0; k<ndata; k++){
      if(k!=n) bin_ave[n] += data[k];
    }
    bin_ave[n] /= ndata - 1.0;
  }

  error = 0.0;
  for(int n=0; n<ndata; n++){ 
    error += (bin_ave[n] - ave) * (bin_ave[n] - ave);
  }
  
  error = sqrt( ndata/(ndata-1.0) * error );

  delete[] bin_ave;
}
