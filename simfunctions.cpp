#include <Rcpp.h>
using namespace Rcpp;


//http://gallery.rcpp.org/articles/stl-random-shuffle/
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
NumericVector randomShuffle(NumericVector a) {
  
  // clone a into b to leave a alone 
  NumericVector b = clone(a);
  
  std::random_shuffle(b.begin(), b.end(), randWrapper);
  
  return b;
}






// [[Rcpp::export]]
void randomTests(){

  //Random number, 0 to 1
  Rcout<<(double)rand()/(double)RAND_MAX<<"\n";
  
}


// [[Rcpp::export]]
double dissimilarityindex(NumericVector pop1, NumericVector pop2){
  
  //Get sums of the two populations
  double sum1 = 0;  
  for(int i = 0; i < pop1.size(); ++i) {
    sum1 += pop1[i];
  }
  
  double sum2 = 0;
  for(int i = 0; i < pop2.size(); ++i) {
    sum2 += pop2[i];
  }
  
  //Get DI
  double DI = 0;  
  
  for(int i = 0; i < pop1.size(); ++i) {
    
    //Absolute...
    //https://github.com/RcppCore/Rcpp/issues/167
    //Works on vectors not doubles unless...
    DI += std::abs((pop1[i]/sum1) - (pop2[i]/sum2));
    
  }
  
  DI *= 0.5;
  
  return DI;
  
}


// [[Rcpp::export]]
List targetDissimilarityIndex(double di_target, double threshold, int breakval,
                              NumericVector pop1, NumericVector pop2, bool keepProportion = false) {
  
  
  //Randomise the population data
  //Get DI each time
  //Break when within threshold distance of target value
  //Or after a set number of iterations
  
  //This will be the final DI as well
  double newDI = 1;
  
  int printDI = 0;
  
  bool withinthreshold = false;
  
  double lastDI = dissimilarityindex(pop1,pop2);
  
  //Get current gap between actual DI and target DI
  double lastgap = std::abs(di_target - lastDI);
  Rcout<<"lastgap at start: "<<lastgap<<"\n";
  
  
  do{
    
    //double lastDI = dissimilarityindex(pop1,pop2);
    
    //Get current gap between actual DI and target DI
    //double lastgap = std::abs(di_target - lastDI);
    
    NumericVector replace_pop1 = clone(pop1);
    NumericVector replace_pop2 = clone(pop2);
  
  
      //Randomise one of the two values
      //Rcout<<rand() % 2;
      if(rand() % 2 == 1) {
        //New random value between 0 and 1
        
        // Rcout<<"1. pop replacement index: "<<rand() % replace_pop1.size()<<"\n";
        // Rcout<<"1. random num: "<<(double)rand() / (double)RAND_MAX<<"\n";
        
        if(keepProportion){
          
          //Maintain population proportion
          //By taking some people from one zone and sticking in another
          //Pick zones to take from and then put in
          //Might sometimes be same zone, doesn't matter
          int swap1 = rand() % replace_pop1.size();
          int swap2 = rand() % replace_pop1.size();
          
          //pull one out
          double swapval1 = replace_pop1[swap1];
          
          //Keep a random proportion of that population
          swapval1 *= (double)rand()/(double)RAND_MAX;
          
          //Subtract from swap1, add to a random zone
          replace_pop1[swap1] -= swapval1;
          replace_pop1[swap2] += swapval1;
          
          //Otherwise we're just randomising, which will be faster
        } else { 
          
          replace_pop1[rand() % replace_pop1.size()] = (double)rand() / (double)RAND_MAX;
          
        }
        
        
        //Work on other population if rand sez so
      } else {
       
        // Rcout<<"2. pop replacement index: "<<rand() % replace_pop1.size()<<"\n";
        // Rcout<<"2. random num: "<<(double)rand() / (double)RAND_MAX<<"\n";
        
        if(keepProportion){
          
          //Maintain population proportion
          //By taking some people from one zone and sticking in another
          //Pick zones to take from and then put in
          //Might sometimes be same zone, doesn't matter
          int swap1 = rand() % replace_pop2.size();
          int swap2 = rand() % replace_pop2.size();
          
          //pull one out
          double swapval1 = replace_pop2[swap1];
          
          //Keep a random proportion of that population
          swapval1 *= (double)rand()/(double)RAND_MAX;
          
          //Subtract from swap1, add to a random zone
          replace_pop2[swap1] -= swapval1;
          replace_pop2[swap2] += swapval1;
          
        } else {
        
          replace_pop2[rand() % replace_pop2.size()] = (double)rand() / (double)RAND_MAX;
          
        }
         
         //End if rand choose 1 of 2 populations 
      }
      
    
    
    newDI = dissimilarityindex(replace_pop1,replace_pop2);
    double newgap = std::abs(di_target - newDI);
    // Rcout<<"lastgap: "<<lastgap<<", newgap: "<<newgap<<"\n";
    
    //Has the random change / population swap made the gap between new DI and target
    //Smaller than gap between original DI and target?
    //If so, keep
    if(newgap < lastgap){
      
      lastgap = newgap;
      
      pop1 = clone(replace_pop1);
      pop2 = clone(replace_pop2);
      
      if(++printDI % 100 == 0){
        Rcout<<newDI<<"\n";
      }
      
    }
    
    if(breakval == 1){
      Rcout<<"Break point reached.\n";
    }
    
    if(newgap < threshold){
      Rcout<<"Threshold reached.\n";
      withinthreshold = true;
    }
    
  }//end do
  while(--breakval > 0 && !withinthreshold);
  
  return List::create(_["pop1"] = pop1,_["pop2"] = pop2);
  
}








// [[Rcpp::export]]
double OLD___getAverageAbsoluteContiguousDifference(NumericVector attribute, int ncol, int nrow){
  
  //For calculating average...
  //Though actually, should be number of neighbours/2 (so not getting A>B and B>A)
  int totalcellnumber = ncol * nrow;
  
  //sum of absolute contig differences
  //Will then divide by zone number to get average ACD
  double totalACD = 0;
  double thisACD;
  
  //Cycle through each attribute value and find neighbouring rook-contiguity cell values
  //Obviously, this only works with a grid cell set-up
  //Also assumes count starts top-left and moves right across cols then down rows
  for(int i = 0; i < attribute.size(); ++i){
    
    //For cell, get absolute contiguous difference on four sides
    
    //LEFT SIDE
    //If on first column, get contiguous value round torus-left
    if(i % ncol == 0){
      thisACD = std::abs(attribute[i] - attribute[i + (ncol-1)]);
      // Rcout<<"i: "<<i<<", left side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + (ncol-1)]<<"\n";
      
    } else {
      thisACD = std::abs(attribute[i] - attribute[i - 1]);
      // Rcout<<"i: "<<i<<", left side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - 1]<<"\n";
      //Rcout<<"i: "<<i<<", left side: "<<thisACD<<"\n";
    }
    
    totalACD += thisACD;
    
    
    //RIGHT SIDE
    //If on last column, get contiguous value round torus-right
    if(i % ncol == ncol-1){
      thisACD = std::abs(attribute[i] - attribute[i - (ncol-1)]);
      // Rcout<<"i: "<<i<<", right side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - (ncol-1)]<<"\n";
      //Rcout<<"i: "<<i<<", right side torus: "<<thisACD<<"\n";
    } else {
      thisACD = std::abs(attribute[i] - attribute[i + 1]);
      // Rcout<<"i: "<<i<<", right side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + 1]<<"\n";
      //Rcout<<"i: "<<i<<", right side: "<<thisACD<<"\n";
    }
    
      
    totalACD += thisACD;
    
    
    //TOP SIDE
    //If on first row, get contiguous value round torus-top
    //If e.g. number of columns is ten and i is between 0 and 9...
    if(i < ncol){
      //eg. 2 + (10 * (10-1)) = 92
      thisACD = std::abs(attribute[i] - attribute[i + ( ncol * (nrow-1) )]);
      // Rcout<<"i: "<<i<<", top side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + ( ncol * (nrow-1) )]<<"\n";
      //Rcout<<"i: "<<i<<", top side torus: "<<thisACD<<"\n";
    } else {
      thisACD = std::abs(attribute[i] - attribute[i - ncol]);
      // Rcout<<"i: "<<i<<", top side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - ncol]<<"\n";
      //Rcout<<"i: "<<i<<", top side: "<<thisACD<<"\n";
    }
    
    totalACD += thisACD;
    
    
    //BOTTOM SIDE
    //If on last row, get contiguous value round torus-bottom
    //if e.g. for 10*10, i > 89 (between 90 and 99)
    if(i > (ncol * (nrow-1))-1 ){
      //eg. 92 - (10 * (10-1)) = 92
      thisACD = std::abs(attribute[i] - attribute[i - ( ncol * (nrow-1) )]);
      // Rcout<<"i: "<<i<<", bottom side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - ( ncol * (nrow-1) )]<<"\n";
      //Rcout<<"i: "<<i<<", bottom side torus: "<<thisACD<<"\n";
    } else {
      thisACD = std::abs(attribute[i] - attribute[i + ncol]);
      // Rcout<<"i: "<<i<<", bottom side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + ncol]<<"\n";
      //Rcout<<"i: "<<i<<", bottom side: "<<thisACD<<"\n";
    }
    
    totalACD += thisACD;
    
    
    
  }
  
  //return per-cell average ACD
  return totalACD / totalcellnumber;
  
}





//Re-writing to allow for an ACD high-pass filter
//(err, only finding values for ACDs above a given cutoff)
// [[Rcpp::export]]
double getAverageAbsoluteContiguousDifference(NumericVector attribute, int ncol, int nrow, double cutoff){
  
  //For calculating average...
  //Though actually, should be number of neighbours/2 (so not getting A>B and B>A)
  //int totalcellnumber = ncol * nrow;
  
  //Rather than totalcellnumber, we can count the actual number of borders
  //And get the average for that
  //As some will be excluded by cutoff, so...
  int totalborders = 0;
  
  //sum of absolute contig differences
  //Will then divide by zone number to get average ACD
  double totalACD = 0;
  double thisACD;
  
  //Cycle through each attribute value and find neighbouring rook-contiguity cell values
  //Obviously, this only works with a grid cell set-up
  //Also assumes count starts top-left and moves right across cols then down rows
  for(int i = 0; i < attribute.size(); ++i){
    
    //For cell, get absolute contiguous difference on four sides
    
    //LEFT SIDE
    //If on first column, get contiguous value round torus-left
    if(i % ncol == 0){
      thisACD = std::abs(attribute[i] - attribute[i + (ncol-1)]);
      // Rcout<<"i: "<<i<<", left side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + (ncol-1)]<<"\n";
      
    } else {
      thisACD = std::abs(attribute[i] - attribute[i - 1]);
      // Rcout<<"i: "<<i<<", left side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - 1]<<"\n";
      //Rcout<<"i: "<<i<<", left side: "<<thisACD<<"\n";
    }
    
    //Keep only if above cutoff
    if(thisACD >= cutoff){
      totalACD += thisACD;
      ++totalborders;
    }
    
    
    //RIGHT SIDE
    //If on last column, get contiguous value round torus-right
    if(i % ncol == ncol-1){
      thisACD = std::abs(attribute[i] - attribute[i - (ncol-1)]);
      // Rcout<<"i: "<<i<<", right side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - (ncol-1)]<<"\n";
      //Rcout<<"i: "<<i<<", right side torus: "<<thisACD<<"\n";
    } else {
      thisACD = std::abs(attribute[i] - attribute[i + 1]);
      // Rcout<<"i: "<<i<<", right side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + 1]<<"\n";
      //Rcout<<"i: "<<i<<", right side: "<<thisACD<<"\n";
    }
    
    
    //Keep only if above cutoff
    if(thisACD >= cutoff){
      totalACD += thisACD;
      ++totalborders;
    }
    
    
    
    //TOP SIDE
    //If on first row, get contiguous value round torus-top
    //If e.g. number of columns is ten and i is between 0 and 9...
    if(i < ncol){
      //eg. 2 + (10 * (10-1)) = 92
      thisACD = std::abs(attribute[i] - attribute[i + ( ncol * (nrow-1) )]);
      // Rcout<<"i: "<<i<<", top side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + ( ncol * (nrow-1) )]<<"\n";
      //Rcout<<"i: "<<i<<", top side torus: "<<thisACD<<"\n";
    } else {
      thisACD = std::abs(attribute[i] - attribute[i - ncol]);
      // Rcout<<"i: "<<i<<", top side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - ncol]<<"\n";
      //Rcout<<"i: "<<i<<", top side: "<<thisACD<<"\n";
    }
    
    //Keep only if above cutoff
    if(thisACD >= cutoff){
      totalACD += thisACD;
      ++totalborders;
    }
    
    
    
    //BOTTOM SIDE
    //If on last row, get contiguous value round torus-bottom
    //if e.g. for 10*10, i > 89 (between 90 and 99)
    if(i > (ncol * (nrow-1))-1 ){
      //eg. 92 - (10 * (10-1)) = 92
      thisACD = std::abs(attribute[i] - attribute[i - ( ncol * (nrow-1) )]);
      // Rcout<<"i: "<<i<<", bottom side torus. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i - ( ncol * (nrow-1) )]<<"\n";
      //Rcout<<"i: "<<i<<", bottom side torus: "<<thisACD<<"\n";
    } else {
      thisACD = std::abs(attribute[i] - attribute[i + ncol]);
      // Rcout<<"i: "<<i<<", bottom side. Attribute here: "<<attribute[i]<<", attribute neighbour: "<<attribute[i + ncol]<<"\n";
      //Rcout<<"i: "<<i<<", bottom side: "<<thisACD<<"\n";
    }
    
    //Keep only if above cutoff
    if(thisACD >= cutoff){
      totalACD += thisACD;
      ++totalborders;
    }
    
    
  }
  
  //return per-cell average ACD
  return totalACD / totalborders;
  
}











//reminder: including secondpop so that 
//Movement of main pop can be matched, repeated for the second population so that DI values remain fixed
// [[Rcpp::export]]
List optimiseAverageAbsoluteContiguousDifference(
    NumericVector attribute, NumericVector secondpop,
    int ncol, int nrow, bool maximise, int breakval, double cutoff){
  
  //bool maximise: if true, maximise otherwise minimise
  
  
  //set up
  double newACD = 1;
  
  int printACD = 0;
  
  double lastACD = getAverageAbsoluteContiguousDifference(attribute,ncol,nrow,cutoff);
  
  
  do{
    
    //copy of attribute to swap zone values around
    NumericVector replace_attr = clone(attribute);
    NumericVector replace_secondpop = clone(secondpop);
    
    //Two indexes for zones to swap values
    //Doesn't matter if we occasionally replace with self
    int swap1 = rand() % replace_attr.size();
    int swap2 = rand() % replace_attr.size();
    
    //pull one out so we can overwrite
    double swapval1 = replace_attr[swap1];
    
    replace_attr[swap1] = replace_attr[swap2];
    replace_attr[swap2] = swapval1;
    
    //Repeat for the second population so that DI values remain fixed
    swapval1 = replace_secondpop[swap1];
    
    replace_secondpop[swap1] = replace_secondpop[swap2];
    replace_secondpop[swap2] = swapval1;
    
    
    
    newACD = getAverageAbsoluteContiguousDifference(replace_attr,ncol,nrow,cutoff);
    
    if(maximise){
      
      if(newACD > lastACD){
        
        lastACD = newACD;
        attribute = clone(replace_attr);
        secondpop = clone(replace_secondpop);
        
        if(++printACD % 100 == 0){
          Rcout<<newACD<<"\n";
        }
        
      }//end if newACD > 
      
    }//end if maximise
    //if minimising
    else {
      
      if(newACD < lastACD){
        
        lastACD = newACD;
        attribute = clone(replace_attr);
        secondpop = clone(replace_secondpop);
        
        if(++printACD % 100 == 0){
          Rcout<<newACD<<"\n";
        }
        
      }//end if newACD > 
      
    }
    
    
    if(breakval == 1){
      Rcout<<"Break point reached.\n";
    }
    
    
    
  }//end do
  while(--breakval > 0);
  //while(--breakval > 0 && !withinthreshold);
  
  //Return both optimised attribute
  //And second pop values tagging along, having been swapped too
  //To keep DI values identical
  
  //return attribute;
  return List::create(_["attribute"] = attribute,_["secondpop"] = secondpop);
  
}






// [[Rcpp::export]]
NumericVector getRepeatedAACDfromPermutedCells(NumericVector attribute, int ncol, int nrow, int numreps, double cutoff){
  
  //Vector to store each results, length of number of reps
  //standard vector may be quicker?
  //https://stackoverflow.com/questions/41602024/should-i-prefer-rcppnumericvector-over-stdvector
  //std::vector<double> results(numreps);
  NumericVector results(numreps);
  
  NumericVector shuffled;
  
  for(int i = 0; i < results.size(); ++i){
    
    shuffled = randomShuffle(attribute);
    
    results[i] = getAverageAbsoluteContiguousDifference(shuffled,ncol,nrow,cutoff);
    
  }
  
  return results; 

}







// [[Rcpp::export]]
void testNeighbours(Rcpp::List nblist, int index){
  
  NumericVector nbs = nblist[index];
  
  for(int i = 0; i < nbs.size(); ++i){
    Rcout<<nbs[i]<<" ";
  }
  
}






//take in cell atttibutes
//And that cell's list of neighbours
// [[Rcpp::export]]
void displayAllNeighbours(NumericVector attribute, Rcpp::List nblist){
  
  //cycle through all attributes
  for(int i = 0; i < attribute.size(); ++i){
  
    NumericVector nbs = nblist[i];
    Rcout<<"i:"<<i<<", neighbour list: ";
  
    //Cycle through all neighbours
    for(int j = 0; j < nbs.size(); ++j){
      Rcout<<nbs[j]<<" ";
    }
    
    Rcout<<"\n";
      
  }
  
}







// [[Rcpp::export]]
std::vector<double> weightedNeighbourIndex(NumericVector attribute, Rcpp::List nblist){
  
  //Vector for storing the actual weighted score
  //Want to look at that
  //https://stackoverflow.com/questions/30129968/how-to-initialize-numericvector-to-a-specific-size-after-it-has-been-declared
  // NumericVector indexscore = NumericVector(attribute.size());
  std::vector<double> indexscore;
  
  
  //cycle through all attributes
  for(int i = 0; i < attribute.size(); ++i){
    
    //Neighbours of i
    NumericVector nbs_of_i = nblist[i];
    
    //Cycle through all of i's neighbours. i and j make dyads/pairs we want to check in turn
    for(int j = 0; j < nbs_of_i.size(); ++j){
      
      // Rcout<<nbs_of_i[j]<<"\n";
      
      //Get neighbours of j (which are neighbours of this particular neighbour of i)
      //SUBTRACT ONE! 
      //The neighbours lists are indexed 1 to n in R.
      //The correct index is going to be nbs_of_i[j]-1 innit?
      NumericVector nbs_of_j = nblist[nbs_of_i[j]-1];
      
      
      
      
      //Add ACD of the found pairs to this
      std::vector<double> acds; 
      
      //cycle through neighbours of neighbours of cell j
      for(int k = 0; k < nbs_of_j.size(); ++k){
        
        NumericVector nbs_of_j_neighbour = nblist[nbs_of_j[k]-1];
        
        //for each of those neighbours of this j neighbour
        //See if it matches any neighbours of i
        //These will be pairs we want to keep
        //(though excluding j itself, which is also a neighbour of neighbour of j)
        for(int m = 0; m < nbs_of_j_neighbour.size(); ++m){
          
          for(int n = 0; n < nbs_of_i.size(); ++n){
            
            if(nbs_of_j_neighbour[m]==nbs_of_i[n]){
              
              //Should exclude i 
              //If we're looking at nbs_of_j when nbs_of_j=i
              //And j when nbs_of_j_neighbour = j
              //as i can be j neighbour and j can be neighbour of j neighbour
              //Or possibly: neither index at m nor n should equal index at i or j. M and N are same value so only need to test one of those.
              //Also have to exclude j neighbour being i
              if(
                ((nbs_of_j_neighbour[m]!=i+1 && nbs_of_j_neighbour[m]!=nbs_of_i[j])&&nbs_of_j[k]!=i+1)
              ){
              
                // Rcout<<"Adjusted to start at 1: i:"<<i+1<<" j:"<<nbs_of_i[j]<<" j_nb:"<<nbs_of_j[k]<<" m: "<<nbs_of_j_neighbour[m]<<" n:"<<nbs_of_i[n]<<"\n";
                
                //Find ACD for that pair and add to vector
                //indexed at  and either nbs_of_j_neighbour[m] or nbs_of_i[n]
                acds.push_back(std::abs(attribute[nbs_of_j[k]-1] - attribute[nbs_of_i[n]-1]));
                
              
              }//end if
              
            }//end if
            
          }//end for n
          
        }//end for m
        
      }//end for k
      
      
      //Add ij ACD to the already-collected ACDs from the (max 2) neighbouring pairs
      //double ij_ACD = std::abs(attribute[i] - attribute[nbs_of_i[j]-1]);
      acds.push_back(std::abs(attribute[i] - attribute[nbs_of_i[j]-1]));
      
      //Add to weighted ACD measure index for this border
      //Find average cos it'll vary in length
      //Total first
      // for(int i = 0; i < acds.size();++i){}
      //https://stackoverflow.com/questions/3221812/how-to-sum-up-elements-of-a-c-vector
      // double sum_of_elems = std::accumulate(acds.begin(), acds.end(), 0);
      
      //That didn't seem to work
      double sum_of_elems = 0;
      //https://stackoverflow.com/questions/7984955/what-is-wrong-with-my-for-loops-i-get-warnings-comparison-between-signed-and-u
      //There are more complications in that thread
      for(unsigned n = 0; n < acds.size(); ++n){
        sum_of_elems += acds[n];
      }
        
      
      // Rcout<<sum_of_elems<<"\n";
      
      indexscore.push_back(sum_of_elems/static_cast<double>(acds.size()));
      
    }//end for j
  }//end for i
  
  return indexscore;
  
}





//Same as above but returning single value, the weighted AACD mean for the whole grid
// [[Rcpp::export]]
double getWeightedNeighbourIndexMean(NumericVector attribute, Rcpp::List nblist){

  //Vector for storing the actual weighted score
  //Want to look at that
  //https://stackoverflow.com/questions/30129968/how-to-initialize-numericvector-to-a-specific-size-after-it-has-been-declared
  // NumericVector indexscore = NumericVector(attribute.size());
  std::vector<double> indexscore;

  //cycle through all attributes
  for(int i = 0; i < attribute.size(); ++i){

    //Neighbours of i
    NumericVector nbs_of_i = nblist[i];

    //Cycle through all of i's neighbours. i and j make dyads/pairs we want to check in turn
    for(int j = 0; j < nbs_of_i.size(); ++j){

      // Rcout<<nbs_of_i[j]<<"\n";

      //Get neighbours of j (which are neighbours of this particular neighbour of i)
      //SUBTRACT ONE!
      //The neighbours lists are indexed 1 to n in R.
      //The correct index is going to be nbs_of_i[j]-1 innit?
      NumericVector nbs_of_j = nblist[nbs_of_i[j]-1];




      //Add ACD of the found pairs to this
      std::vector<double> acds;

      //cycle through neighbours of neighbours of cell j
      for(int k = 0; k < nbs_of_j.size(); ++k){

        NumericVector nbs_of_j_neighbour = nblist[nbs_of_j[k]-1];

        //for each of those neighbours of this j neighbour
        //See if it matches any neighbours of i
        //These will be pairs we want to keep
        //(though excluding j itself, which is also a neighbour of neighbour of j)
        for(int m = 0; m < nbs_of_j_neighbour.size(); ++m){

          for(int n = 0; n < nbs_of_i.size(); ++n){

            if(nbs_of_j_neighbour[m]==nbs_of_i[n]){

              //Should exclude i
              //If we're looking at nbs_of_j when nbs_of_j=i
              //And j when nbs_of_j_neighbour = j
              //as i can be j neighbour and j can be neighbour of j neighbour
              //Or possibly: neither index at m nor n should equal index at i or j. M and N are same value so only need to test one of those.
              //Also have to exclude j neighbour being i
              if(
                ((nbs_of_j_neighbour[m]!=i+1 && nbs_of_j_neighbour[m]!=nbs_of_i[j])&&nbs_of_j[k]!=i+1)
              ){

                // Rcout<<"Adjusted to start at 1: i:"<<i+1<<" j:"<<nbs_of_i[j]<<" j_nb:"<<nbs_of_j[k]<<" m: "<<nbs_of_j_neighbour[m]<<" n:"<<nbs_of_i[n]<<"\n";

                //Find ACD for that pair and add to vector
                //indexed at  and either nbs_of_j_neighbour[m] or nbs_of_i[n]
                acds.push_back(std::abs(attribute[nbs_of_j[k]-1] - attribute[nbs_of_i[n]-1]));


              }//end if

            }//end if

          }//end for n

        }//end for m

      }//end for k


      //Add ij ACD to the already-collected ACDs from the (max 2) neighbouring pairs
      //double ij_ACD = std::abs(attribute[i] - attribute[nbs_of_i[j]-1]);
      acds.push_back(std::abs(attribute[i] - attribute[nbs_of_i[j]-1]));

      //Add to weighted ACD measure index for this border
      //Find average cos it'll vary in length
      //Total first
      // for(int i = 0; i < acds.size();++i){}
      //https://stackoverflow.com/questions/3221812/how-to-sum-up-elements-of-a-c-vector
      // double sum_of_elems = std::accumulate(acds.begin(), acds.end(), 0);

      //That didn't seem to work
      double sum_of_elems = 0;
      //https://stackoverflow.com/questions/7984955/what-is-wrong-with-my-for-loops-i-get-warnings-comparison-between-signed-and-u
      //There are more complications in that thread
      for(unsigned n = 0; n < acds.size(); ++n){
        sum_of_elems += acds[n];
      }


      // Rcout<<sum_of_elems<<"\n";

      indexscore.push_back(sum_of_elems/static_cast<double>(acds.size()));

    }//end for j
  }//end for i


  //get mean index score
  double sum_of_elems = 0;
  for(unsigned n = 0; n < indexscore.size(); ++n){
    sum_of_elems += indexscore[n];
  }

  return sum_of_elems / static_cast<double>(indexscore.size());

}












//Hack until we make these generic
// [[Rcpp::export]]
List optimiseWEIGHTEDAverageAbsoluteContiguousDifference(
    NumericVector attribute, NumericVector secondpop, Rcpp::List nblist,
     bool maximise, int breakval, double cutoff){

  //bool maximise: if true, maximise otherwise minimise


  //set up
  double newACD = 1;

  int printACD = 0;

  double lastACD = getWeightedNeighbourIndexMean(attribute,nblist);


  do{

    //copy of attribute to swap zone values around
    NumericVector replace_attr = clone(attribute);
    NumericVector replace_secondpop = clone(secondpop);

    //Two indexes for zones to swap values
    //Doesn't matter if we occasionally replace with self
    int swap1 = rand() % replace_attr.size();
    int swap2 = rand() % replace_attr.size();

    //pull one out so we can overwrite
    double swapval1 = replace_attr[swap1];

    replace_attr[swap1] = replace_attr[swap2];
    replace_attr[swap2] = swapval1;

    //Repeat for the second population so that DI values remain fixed
    swapval1 = replace_secondpop[swap1];

    replace_secondpop[swap1] = replace_secondpop[swap2];
    replace_secondpop[swap2] = swapval1;



    newACD = getWeightedNeighbourIndexMean(replace_attr,nblist);

    if(maximise){

      if(newACD > lastACD){

        lastACD = newACD;
        attribute = clone(replace_attr);
        secondpop = clone(replace_secondpop);

        if(++printACD % 200 == 0){
          Rcout<<breakval<<": "<<newACD<<"\n";
        }

      }//end if newACD >

    }//end if maximise
    //if minimising
    else {

      if(newACD < lastACD){

        lastACD = newACD;
        attribute = clone(replace_attr);
        secondpop = clone(replace_secondpop);

        if(++printACD % 100 == 0){
          Rcout<<breakval<<": "<<newACD<<"\n";
        }

      }//end if newACD >

    }


    if(breakval == 1){
      Rcout<<"Break point reached.\n";
    }



  }//end do
  while(--breakval > 0);
  //while(--breakval > 0 && !withinthreshold);

  //Return both optimised attribute
  //And second pop values tagging along, having been swapped too
  //To keep DI values identical

  //return attribute;
  return List::create(_["attribute"] = attribute,_["secondpop"] = secondpop);

}


