#pragma once
#include<vector>
#include<algorithm>

class BS
{
public:
	BS();
	~BS();

	void initialize(int n);

	int sum(void);

private:

	std::vector<int> arr;

};


