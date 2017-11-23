//***********************************************************************************************
//              Project 1 : Indrect branch predictor Implementation using VPC 
//              
//              Conditional Branch predictor used : Gshare 
//
//              By : Srinivas Kovuri  UIN:826009011
//                   Texas A & M University, email - srinivaskovuri@tamu.edu 
//***********************************************************************************************
// my_predictor.h
// This file contains a sample my_predictor class.
// It has a simple 32,768-entry gshare with a history length of 15 and a
// simple direct-mapped branch target buffer for indirect branch prediction.
static unsigned int vpc_gen[32] = { 65280,61680,52428,4390, 30735, 255, 3855, 4080, 2735, 6253, 56305, 24864, 52280, 42300, 46198, 2212, 13739, 5471, 9941, 15962,55943,22507,32231,44979,2630,14315,27870,35864,14301,31866,51083,3260};    //1111111100000000 , 1111000011110000 , 1100110011001100 , 1010101010101010, 1111000000001111, 0000000011111111, 0000111100001111, 1111 1111 0000, 101010101111 ...  //made an educated guess for choosing intial values to avoid aliasing but as the no. of entries increased it has become difficult to come up with the values. Used random generator for later values 

class my_update : public branch_update {
public:
	unsigned long int index;
};

class my_predictor : public branch_predictor {
public:
	#define HISTORY_LENGTH	16
	#define TABLE_BITS	16
	#define NUM_ITER         32     // Max number iterations to predict the target or train the predictor in worst case.
 
        //For Piecewise Linear Conditional Branch Predictor
	#define GA_LENGTH       8       //If this is 1, it will be preceptron based predictor
        #define N_LENGTH        8       //If this is 1, it will be path based predictor
	int output;
        bool predict_output;
        int W[1<<N_LENGTH][(1<<GA_LENGTH)][HISTORY_LENGTH];
        unsigned char GA[(HISTORY_LENGTH)];

	my_update u;
	branch_info bi;
	unsigned int history;
	unsigned char tab[1<<TABLE_BITS];
	//unsigned long int BTB[1<<(TABLE_BITS+NUM_ITER/3)];  //array structure as BTB. Index corresponds to virtual PC and the value to its corresponding target
	unsigned long int BTB[1<<TABLE_BITS];  //array structure as BTB. Index corresponds to virtual PC and the value to its corresponding target
    	unsigned int vpc;  //virtual PC 
   	unsigned int vghr;   //virtual GHR
    	int predict_iter;   //no of iterations took for predicting the target which will be used for training VPC
    	unsigned int predict_target;   //for storing the predicted target address
	my_predictor (void) : history(0),vpc(0),vghr(0),predict_target(0) { 
		memset (tab, 0, sizeof (tab));
        	//set add initialisations here
	        memset (BTB,0, sizeof(BTB)); 
		memset (W,0,sizeof(W));
                memset (GA,0,sizeof(GA)); 
	}
	
	void update_conditional_LPCBP(unsigned int address, bool bTaken){     //for training Linear Piecewise Conditional Predictor
                unsigned int Theta = 2.14*(HISTORY_LENGTH)+20.8;
                address = address & ((1<<N_LENGTH)-1);
                predict_output = u.direction_prediction();
                if((output < Theta && output > (-1*Theta) ) || predict_output != bTaken){
                        if(bTaken){
                                if(W[address][0][0] < 127){
                                        W[address][0][0] += 1;
                                }
                        }
                        else{
                                if(W[address][0][0] > -127){
                                        W[address][0][0] -= 1;
                                }
                        }
                        for(int i=1; i < (HISTORY_LENGTH) + 1; i++){
                                if(((history >> i) & 1)){
                                        if(W[address][GA[i]][i] < 127){
                                                W[address][GA[i]%(1<<GA_LENGTH)][i] += 1;
                                        }
                                }
                                else{
                                        if(W[address][GA[i]][i] > -127){
                                                W[address][GA[i]%(1<<GA_LENGTH)][i] -= 1;
                                        }
                                }
                        }
                }
                for(int i=HISTORY_LENGTH; i >1; i--){
                        GA[i] = (GA[i-1]) % (1<< GA_LENGTH);  //shifting the address to accomodate the latest entry in GA to indice 1 and as you move up the indice you get recently traversed addresses
                }
                GA[1] = address & ((1 << GA_LENGTH)-1);
                history <<= 1;
                history |= bTaken;
                history &= ((1 << HISTORY_LENGTH)-1);    //logic for updating the Global history register

        }

	void predict_conditional_LPCBP(unsigned int address){   //for predicting Linear piecewise conditional Branch predictor
                address = address & ((1<<N_LENGTH)-1);  //masking address to fall within the range of m defined for piecewise linear predictor
                output = W[address][0][0];
                for(int i=1; i <= HISTORY_LENGTH; i++){
                        if((history >> (i-1)) & 1){  //if the bit within the GHR is taken update output with adding weight else subtract it based on saturation
                                if(output < 127){
                                        output +=  W[address][GA[i]%(1<<GA_LENGTH)][i];
                                }
                        }
                        else{
                                if(output > -127){
                                        output -=  W[address][GA[i]%(1<<GA_LENGTH)][i];
                                }
                        }
                }
                if(output >= 0){
                        u.direction_prediction(true);   //if the output >= 0 predict the direction outcome as taken
                }
                else{
                        u.direction_prediction(false);  //predicting the direction outcome as false
                }
        }

	void update_conditional_BP(unsigned int virtualPC, unsigned int virtualGHR, bool bTaken){    //procedure for updating the GShare conditional branch predictor  used within the training algorithm
		unsigned long int index = (virtualGHR ^ (virtualPC & ((1 << TABLE_BITS)-1))) & ((1 << TABLE_BITS)-1);
		history <<= 1;
		history |= bTaken;
		history &= (1<<HISTORY_LENGTH)-1;
		if(bTaken){
			if(tab[index] <3 ){
				tab[index]++;
			}
		}
		else{
			if(tab[index] > 0 ){
				tab[index]--;
			}
		}		
	}

	branch_update *predict (branch_info & b) {
        	bi = b;                
		if (b.br_flags & BR_CONDITIONAL) {
			//predict_conditional_LPCBP(b.address);  //for predicting Linear Piecewise Conditional branch predictor
			
			//For predicting GShare 
			u.index = 
				  ((history << (TABLE_BITS - HISTORY_LENGTH)) 
				^ (b.address & ((1<<TABLE_BITS)-1)))&((1<<TABLE_BITS)-1);      // calculating the index by XORing the GHR with PC
			u.direction_prediction (tab[u.index] >> 1);     //predicting the direction. 3,2 -> taken & 1,0 not taken
		} else {
			u.direction_prediction (true);
		}
		if (b.br_flags & BR_INDIRECT) {
		//	Implementing the VPC algorithm code here for predicting Indirect branch target
			int iter=0;
            		vpc  = b.address & ((1<<TABLE_BITS)-1);   //initialising vpc <- Indirect branch PC
			
			// Disable this for using Updates to VPC training algorithm using GShare and Enable when using Piecewise Linear Predictor
			if(0){
				for(iter=0;iter<NUM_ITER;iter++){
					predict_iter = iter;
                                	if(BTB[vpc] != 0 ){   //Continue iterating the VPC predicting algorithm only there is a BTB hit
                                        	predict_conditional_LPCBP(vpc);
                                        	if(u.direction_prediction()){
                                                	//to point next PC to the branch target if there is a BTB hit and the conditional branch is taken
                                                	u.target_prediction(BTB[vpc]);
                                               		predict_target = u.target_prediction();
							break; //stop iterating when the Conditional branch is taken and BTB hit
						}
					}
					else{
			        		u.target_prediction(false);  
						break;     //stop iterating as there is no BTB hit
					}
					vpc = (b.address ^ (vpc_gen[iter])) & ((1 << TABLE_BITS) - 1);   //Generating different VPCs by XORing with hard-coded values with PC
				}
				if(iter == NUM_ITER){  
					u.target_prediction(false);  //when there is no "TAKEN" outcome for conditional branch
				}
			}
			
			// Enable for using Updates to VPC training algorithm using GShare
			if(1){
				vghr = history;                           // vpc <- GHR
            			for(iter=0; iter < NUM_ITER ; iter++){
		        		predict_iter = iter;
		        		if(BTB[vpc] != 0 ){   //Continue iterating the VPC predicting algorithm only there is a BTB hit
						u.index = (vpc ^ vghr) & ((1 << TABLE_BITS)-1);
						u.direction_prediction(tab[u.index] >> 1);	
						if(u.direction_prediction()){    //if the outcome of vpc is taken 
     							u.target_prediction(BTB[vpc]);
							predict_target = u.target_prediction();
							break;  //stop iterating when the Conditional branch is taken and BTB hit for a corresponding vghr and vpc
						}
					}
				        else{
				        	u.target_prediction(false);  
						break;     //stop iterating as there is no BTB hit
				        }
					vpc = (b.address ^ (vpc_gen[iter])) & ((1 << TABLE_BITS) - 1);   //Generating different VPCs by XORing with hard-coded values with PC
					vghr = (vghr << 1) & ((1 << TABLE_BITS) - 1);   //left shifting VGHR
				}
				
				if(iter == NUM_ITER){  
					u.target_prediction(false);  //when there is no "TAKEN" outcome for conditional branch
				}
			}
		}
		return &u;
	}

	void update (branch_update *u, bool taken, unsigned int target) {
		unsigned int vpc_BTB_Miss = 0;
		unsigned int vghr_BTB_Miss = 0;
		if (bi.br_flags & BR_CONDITIONAL) {
			unsigned char *c = &tab[((my_update*)u)->index];
			if (taken) {
				if (*c < 3) (*c)++;
			} else {
				if (*c > 0) (*c)--;
			}
			history <<= 1;
			history |= taken;
			history &= (1<<HISTORY_LENGTH)-1;
		}
		if (bi.br_flags & BR_INDIRECT) {
			//Implementing the VPC algorithm for training
			// 1 . To update when VPC target prediction is correct
            		vpc  = bi.address & ((1 << TABLE_BITS)-1);
			vghr = history;
			if(BTB[vpc] != 0){           //When there is a BTB hit for a vpc
				if(predict_target == target) {   //if the predicted target is correct train the BTB and Conditional predictor
					for (int iter=0; iter < predict_iter+1 ; iter ++){
						if(iter == predict_iter){
							//update_conditional_LPCBP(vpc,true);    //update the condition BP with taken for correctly predicted vpca
							update_conditional_BP(vpc,vghr,true);    //update the condition BP with taken for correctly predicted vpca
						}
						else{
							//update_conditional_LPCBP(vpc,false);    //update the condition BP when not taken 
							update_conditional_BP(vpc,vghr,false);   ////update the condition BP with not taken
						}
						vpc = (bi.address ^ vpc_gen[iter]) & ((1 << TABLE_BITS) - 1);  //Generating different VPCs by XORing with hard-coded values with PC
						vghr = (vghr << 1) & ((1 << TABLE_BITS) - 1);  //VGHR <- VGHR <<1 
					}
				}
				else{   //when target address is mispredicted
					int iter;
                    			for (iter=0 ; iter < NUM_ITER ; iter++){   //To iterate till the correct target address found and update corresponding conditional predictor with taken/Not taken accordingly (1. For correct target address update as taken and not taken for rest that are traversed )
						if(BTB[vpc] != 0 ){    // If there is a BTB Hit
							u->target_prediction(BTB[vpc]);
							predict_target = u->target_prediction();
							if(predict_target == target){
								//update_conditional_LPCBP(vpc,true);    //update the condition BP when not taken 
								update_conditional_BP(vpc,vghr,true);   
								break;
							}
							else {
								//update_conditional_LPCBP(vpc,false);    //update the condition BP when not taken 
								update_conditional_BP(vpc,vghr,false);								
							}
			        			vpc = ((bi.address & ((1<<TABLE_BITS)-1)) ^ (vpc_gen[iter] & ((1 << TABLE_BITS)-1))); 
			        			vghr = (vghr << 1) & ((1 << TABLE_BITS)-1);  
						}
						else{   //store the last vpca & vghr for last BTB Miss to update it as correct target  incase no  target found and update it as taken for Conditional BP
							vpc_BTB_Miss = vpc;    
							vghr_BTB_Miss = vghr;
						}
					}
					if(iter == NUM_ITER){
						BTB[vpc_BTB_Miss] = target;      //Insert correct target within the BTB for the last BTB Miss occured
						update_conditional_BP(vpc_BTB_Miss,vghr_BTB_Miss,true);		//update the direction of Conditional branch of virtual PC with BTB Miss as true		
						//update_conditional_LPCBP(vpc_BTB_Miss,true);    //update the condition BP when not taken 
					}
				}
			}
			else{    //when there is a BTB miss for the first PC, insert the target for a virtual PC within BTB and update as taken within Conditional BP
				BTB[vpc] = target;
				update_conditional_BP(vpc,vghr,true);			
				//update_conditional_LPCBP(vpc,true);			
			}
		}
	}
};
