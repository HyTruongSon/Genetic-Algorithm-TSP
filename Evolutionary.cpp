// Software: Genetic (Evolutionary) Algorithm to solve Travelling Saleman Problem (TSP)
// Author: Hy Truong Son
// Major: BSc. Computer Science
// Class: 2013 - 2016
// Institution: Eotvos Lorand University
// Email: sonpascal93@gmail.com
// Website: http://people.inf.elte.hu/hytruongson/
// Final update: October 4th, 2015
// Copyright 2015 (c) Hy Truong Son. All rights reserved. Only use for academic purposes.

#include <iostream>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <iterator>
#include <algorithm>
#include <fstream>

using namespace std;

const int MaxN = 2000;
const int Max_nPopulation = 20000;

int nIterations = 100;
int nPopulation = 1000;
int Selection_Percent = 30;
int Mutation_nSwaps = 10;

string FileName;
int nPoints;
double dist[MaxN + 1][MaxN + 1];

struct aPoint {
	int x, y;
	
	double FindLen(const aPoint another){
		int dx = (another.x - x) * (another.x - x);
		int dy = (another.y - y) * (another.y - y);
		return sqrt(dx + dy);
	}
};
aPoint Point[MaxN + 1];

int temp[MaxN];
bool exist[MaxN + 1];

struct aGene {
	double total;
	int info[MaxN];
	
	bool operator < (const aGene another) const{
		if (total < another.total) return true;
		return false;
	}
	
	void calculate_total(){
		total = 0.0;
		for (int i = 0; i < nPoints; i++){
			int u = info[i];
			int v = info[(i + 1) % nPoints];
			total += dist[u][v];
		}
	}
	
	void random(){
		for (int i = 0; i < nPoints; i++) info[i] = i + 1;
		
		for (int i = 1; i <= 10 * nPoints; i++){
			int u = rand() % nPoints;
			int v = rand() % nPoints;
			swap(info[u], info[v]);
		}
		
		calculate_total();
	}
	
	void redirect(int i, int j){
		int a = info[i];
		int b = info[(i + 1) % nPoints];
		int c = info[j];
		int d = info[(j + 1) % nPoints];
		
		int count = 0;
		
		int v = j;
		while (true){
			temp[count] = info[v];
			count++;
			if (info[v] == b) break;
			v = (v - 1 + nPoints) % nPoints;
		}
		
		v = (j + 1) % nPoints;
		while (true){
			temp[count] = info[v];
			count++;
			if (info[v] == a) break;
			v = (v + 1) % nPoints;
		}
		
		for (int i = 0; i < nPoints; i++) info[i] = temp[i];
	}
	
	void local_optimization(){
		bool stop = false;
		
		while (true){
			stop = true;
			
			for (int i = 0; i < nPoints; i++){
				int a = info[i];
				int b = info[(i + 1) % nPoints];
			
				int j = (i + 2) % nPoints;
				while (true){
					int c = info[j];
					int d = info[(j + 1) % nPoints];
				
					if (d == a) break;
					
					if (dist[a][c] + dist[b][d] < dist[a][b] + dist[c][d]){
						stop = false;
						redirect(i, j);
						break;
					}
				
					j = (j + 1) % nPoints;
				}
				
				if (!stop) break;
			}
			
			if (stop) break;
		}
		
		calculate_total();
	}
	
	void mutation(){
		for (int i = 0; i < Mutation_nSwaps; i++){
			int u = rand() % nPoints;
			int v = rand() % nPoints;
			swap(info[u], info[v]);
		}
		
		local_optimization();
	}
};

aGene Gene[Max_nPopulation + 1];
aGene Child[Max_nPopulation + 1];
aGene next[Max_nPopulation + 1];

void Input(){
	string s, temp, s1, s2, s3;
	ifstream file(FileName.c_str(), ios::in);
	
	while (true){
		getline(file, s);
		if (s[0] == 'E') break;
		
		int i = 0;
		int count = 0;
		
		while (i < s.length()){
			if (s[i] == ' '){
				i++;
				continue;
			}
			
			temp = "";
			for (int j = i; j < s.length(); j++)
				if (s[j] != ' '){
					i = j;
					temp += s[j];
				}else break;
			
			count++;
			if (count == 1) s1 = temp; else
				if (count == 2) s2 = temp; else
					s3 = temp;
			i++;
		}
		
		nPoints++;
		i = atoi(s1.c_str());
		Point[i].x = atoi(s2.c_str());
		Point[i].y = atoi(s3.c_str());
	}
	
	file.close();
}

void Init(){
	for (int i = 1; i < nPoints; i++)
		for (int j = i + 1; j <= nPoints; j++){
			dist[i][j] = Point[i].FindLen(Point[j]);
			dist[j][i] = dist[i][j];
		}
		
	for (int i = 0; i < nPopulation; i++){
		Gene[i].random();
		Gene[i].local_optimization();
	}
}

aGene crossover(const aGene Gene1, const aGene Gene2){
	int i = rand() % nPoints;
	
	int j = -1;
	for (int t = 0; t < nPoints; t++)
		if (Gene2.info[t] == Gene1.info[i]){
			j = t;
			break;
		}
	
	memset(exist, false, sizeof(exist));
	aGene res;
	res.info[i] = Gene1.info[i];
	exist[res.info[i]] = true;
	
	int u = (i - 1 + nPoints) % nPoints;
	int v = (i + 1) % nPoints;
	
	i = (i - 1 + nPoints) % nPoints;
	j = (j + 1) % nPoints;
	
	for (int t = 1; t < nPoints; t++)
		if (t % 2 == 0){
			while (true){
				if (exist[Gene1.info[i]]){
					i = (i - 1 + nPoints) % nPoints;
					continue;
				}
				
				res.info[u] = Gene1.info[i];
				exist[Gene1.info[i]] = true;
				
				u = (u - 1 + nPoints) % nPoints;
				i = (i - 1 + nPoints) % nPoints;
				break;
			}
		}else
			while (true){
				if (exist[Gene2.info[j]]){
					j = (j + 1) % nPoints;
					continue;
				}
				
				res.info[v] = Gene2.info[j];
				exist[Gene2.info[j]] = true;
				
				v = (v + 1) % nPoints;
				j = (j + 1) % nPoints;
				break;
			}
	
	res.local_optimization();
	return res;
}

void output(aGene BestCurrent){
	ofstream file("solution.dat");
	file << "Minimum distance found by Evolutionary algorithm: " << BestCurrent.total << endl;
	for (int i = 0; i < nPoints; i++)
		file << BestCurrent.info[i] << endl;
	file.close();
}

void Evolutionary_Algorithm(){
	sort(Gene, Gene + nPopulation);
	
	cout << endl << "Iteration 0" << endl;
	cout << "Best gene: " << Gene[0].total << endl;
	
	for (int iter = 1; iter <= nIterations; iter++){
		cout << endl << "Iteration " << iter << endl;
		
		int nChilds = 0;
		for (int i = 0; i < nPopulation - 1; i++){
			Child[nChilds] = crossover(Gene[i], Gene[i + 1]);
			//Child[nChilds].mutation();
			nChilds++;
			
			Child[nChilds] = crossover(Gene[i + 1], Gene[i]);
			//Child[nChilds].mutation();
			nChilds++;
			
			if ((i + 1) * 100 >= 0.5 * Selection_Percent * nPopulation) break;
		}
		
		for (int i = nPopulation - 1; i > 0; i--){
			Child[nChilds] = crossover(Gene[i], Gene[i - 1]);
			//Child[nChilds].mutation();
			nChilds++;
			
			Child[nChilds] = crossover(Gene[i - 1], Gene[i]);
			//Child[nChilds].mutation();
			nChilds++;
			
			if ((nPopulation - i) * 100 >= 0.5 * Selection_Percent * nPopulation) break;
		}
		
		sort(Child, Child + nChilds);
		
		int i = 0;
		int j = 0;
		
		for (int v = 0; v < nPopulation; v++){
			if (j == nChilds){
				next[v] = Gene[i];
				i++;
				continue;
			}
			
			if (i == nPopulation){
				next[v] = Child[j];
				j++;
				continue;
			}
			
			if (Gene[i] < Child[j]){
				next[v] = Gene[i];
				i++;
			}else{
				next[v] = Child[j];
				j++;
			}
		}
		
		for (int v = 0; v < nPopulation; v++) Gene[v] = next[v];
		cout << "Best gene: " << Gene[0].total << endl;
		
		if (iter % 10 == 0){
			output(Gene[0]);
			cout << endl << "Saved the best current solution to file solution.dat" << endl;
		}
	}
}

void Getting_Parameters(){
	cout << "Data file: ";
	getline(cin, FileName);
	
	string s;
	
	cout << "Number of iterations (default " << nIterations << "): ";
	getline(cin, s);
	if (s.length() > 0)
		nIterations = atoi(s.c_str());
		
	cout << "Population size (default " << nPopulation << "): ";
	getline(cin, s);
	if (s.length() > 0)
		nPopulation = atoi(s.c_str());
		
	cout << "Selection percent (default " << Selection_Percent << "): ";
	getline(cin, s);
	if (s.length() > 0)
		Selection_Percent = atoi(s.c_str());
}

int main(){
	Getting_Parameters();
	Input();
	Init();
	Evolutionary_Algorithm();
	
	output(Gene[0]);
	cout << endl << "Best solution found by Evolutionary Algorithm: " << Gene[0].total << endl;
}
