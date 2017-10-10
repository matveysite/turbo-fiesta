#pragma once
#include <iostream>

class Set {
public:
	Set(int);
	void insert(Set&, int);
	bool check(int) const;
	void add(int);
	void remove(int);
	int getCount() const;
	int getN() const;
	void inverse();
	friend std::ostream& operator<<(std::ostream&, Set&);
	Set& operator=(Set&);
	~Set();
private:
	int n, count;
	bool* set;
};