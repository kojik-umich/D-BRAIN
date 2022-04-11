#include "BS.h"

BS::BS()
{
}

BS::~BS()
{
}

void BS::initialize(int n)
{
	this->arr.resize(n);

	for (int i = 0; i < n; i++)
		this->arr[i] = i;

	return;
}

int BS::sum(void)
{
	int s = 0;

	for (int& i : this->arr)
		s += i;

	return s;
}
