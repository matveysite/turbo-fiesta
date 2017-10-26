#include "Set.h"

Set::Set(int n):count(0), n(n), set(new bool[n]()){}

void Set::insert(Set& set, int k) {
	for(int i = 0; (i < n) && (i < set.n) && (i < k); i++) {
		this->set[i] = set.set[i];
	}
}

Set& Set::operator=(Set& set){
	if(n!=set.n) {
		delete[] this->set;
		n = set.n;
		this->set = new bool[n];
	}
	for(int i = 0; i < n; i++) {
		this->set[i] = set.set[i];
	}
	return *this;
}

bool Set::check(int i) const {
	return set[i];
}

void Set::add(int i) {
	if (!set[i]) {
		count++;
	}
	set[i] = true;
}

void Set::remove(int i) {
	if (set[i]) {
		count--;
	}
	set[i] = false;
}

int Set::getCount() const {
	return count;
}

int Set::getN() const {
	return n;
}

void Set::inverse() {
	count = n - count;
	for (int i = 0; i < n; i++) {
		set[i]= !set[i];
	}
}

Set::~Set() {
	delete[] set;
}

std::ostream & operator<<(std::ostream& out, Set& set) {
	for (int i = 0; i < set.n; i++) {
		if (set.set[i])
			std::wcout << i + 1 << " ";
	}
	return out;
}
