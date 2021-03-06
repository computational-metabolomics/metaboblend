/*
 * MatchListener.h
 *
 *  Created on: Aug 3, 2012
 *      Author: vbonnici
 */
/*
Copyright (c) 2014 by Rosalba Giugno

This library contains portions of other open source products covered by separate
licenses. Please see the corresponding source files for specific terms.

RI is provided under the terms of The MIT License (MIT):

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef MATCHLISTENER_H_
#define MATCHLISTENER_H_

#include <set>
#include <vector>
#include <stdlib.h>

namespace rilib{

class MatchListener{
public:
	long matchcount;
	MatchListener(){
		matchcount = 0;
	}
	virtual ~MatchListener(){};
	virtual void match(int n, int* qIDs, int* rIDs)=0;
};


class EmptyMatchListener : public MatchListener {
public:
	EmptyMatchListener() : MatchListener(){
	}
	virtual void match(int n, int* qIDs, int* rIDs){
		matchcount++;
	};
};

class ConsoleMatchListener : public MatchListener {
public:
	ConsoleMatchListener() : MatchListener(){
	}
	virtual void match(int n, int* qIDs, int* rIDs){
		matchcount++;
		std::cout<< "{";
		for(int i=0; i<n; i++){
			std::cout<< "("<< qIDs[i] <<","<< rIDs[i] <<")";
		}
		std::cout<< "}\n";
	}
};



class NonRedundantMatchListener : public MatchListener {
public:


	std::set< std::vector<int> > matches;

	int j;

	NonRedundantMatchListener() : MatchListener(){
		j = -1;
	}
	virtual void match(int n, int* qIDs, int* rIDs){
		//match_t m(n,qIDs, rIDs);
		//matches.insert(m);

		std::vector<int> v(n);
		for(int i=0; i<n; i++){
			v[ qIDs[i] ] = rIDs[i];
		}

		if(matches.insert(v).second) {
			
			matchcount++;

			std::cout<<"[";
			for(int i=0; i<n-1; i++){
				std::cout << v[i] << ",";
			}
			std::cout<<v[n-1];
			std::cout<<"]\n";
		}
	
		if ( v[0] != j ){
			//std::cout << "starting integer " << v[0] << std::endl;
			//std::cout << "size (start) " << matches.size() << std::endl;
			//std::cout << v[0] << " -- " << j << std::endl;
			matches.clear();
			//if (j != -1){
			//	exit (EXIT_FAILURE);
			//}
			j = v[0];

		}
	}
};


}


#endif /* MATCHLISTENER_H_ */
