% estimation g,g+2,g+4 cycles and
% connected cycles to variables and check nodes i
function [g_ Ng_ Ng2_ Ng4_ Ng_per_u_ Ng2_per_u_ Ng4_per_u_ Ng_per_p_ Ng2_per_p_ Ng4_per_p_ ] = girth(G) 

g_ = 0; 
Ng_  = 0; 
Ng2_ = 0; 
Ng4_ = 0; 

%define the edge matrix from the given graph G 
E = sparse(G); 
Et = sparse(E'); 
 
%let W_len =|W| cardinality of one set 
W_len = length(G(1,:)); 
Ng_per_p_=zeros(size(G(1,:)));   % variable node in cycle
Ng2_per_p_=zeros(size(G(1,:)));  
Ng4_per_p_=zeros(size(G(1,:))); 
%let U_len =|U| cardinality of the other set 
U_len = length(G(:,1)); 
Ng_per_u_=zeros(size(G(:,1))); 
Ng2_per_u_=zeros(size(G(:,1)));   % check node in cycle
Ng4_per_u_=zeros(size(G(:,1))); 

%check |U| and |W|  
if (U_len > W_len) 
    E = Et; 
    Et = G; 
    temp = U_len; 
    U_len = W_len; 
    W_len = temp; 
end 
 
%Compute the recursion base for U Matrix 
P_U_2_      = (E*Et)-diag(diag(E*Et)); 
P_U_2_c2_   = ((P_U_2_).*((P_U_2_)-1)./2); 
 
L_U_0_2_m1_ = max((diag(diag(E*Et))-1),0); 
L_U_0_2_m2_ = max((diag(diag(E*Et))-2),0); 
Iu          = diag(diag(ones(U_len))); 
 
%Compute the recursion base for W Matrix 
P_W_2_      = (Et*E)-diag(diag(Et*E)); 
P_W_2_c2_   = ((P_W_2_).*((P_W_2_)-1)./2); 
 
L_W_0_2_m1_ = max((diag(diag(Et*E))-1),0); 
L_W_0_2_m2_ = max((diag(diag(Et*E))-2),0); 
Iw          = diag(diag(ones(W_len))); 
 
%Count 4-Cycles------------------------------------------------------------ 
 
    %calculate --- L(1,2) 
    L_U_1_2_    = E*L_W_0_2_m1_; 
    L_W_1_2_    = Et*L_U_0_2_m1_; 
 
    %calculate --- P(g-1) 
    P_U_3_      = P_U_2_*E - L_U_1_2_; 
    P_W_3_      = P_U_3_'; 
 
    %calculate --- L(0,g) 
    L_U_0_4_	= P_U_3_*Et.*Iu; 
    L_W_0_4_    = P_W_3_*E.*Iw; 
 
    if (trace(L_U_0_4_) ~= 0) 
        g_  = 4; 
        Ng_ = (trace(L_U_0_4_) / 4); 
        %connection of cycles 4 to every variable nodes
        Ng_per_p_=diag(L_W_0_4_/4);
        %connection of cycles 4 to every check nodes
        Ng_per_u_=diag(L_U_0_4_/4);
  
                
    end 
 
%Count 6-Cycles------------------------------------------------------------ 
 
    %calculate --- L(2,2) 
    L_U_2_2_ = E*L_W_1_2_- diag(diag(E*L_W_1_2_)); 
    L_W_2_2_ = Et*L_U_1_2_- diag(diag(Et*L_U_1_2_)); 
 
    %calculate --- P(g) 
    P_U_4_ = P_U_3_*Et - L_U_0_4_ - L_U_2_2_; 
    P_W_4_ = P_W_3_*E - L_W_0_4_ - L_W_2_2_; 
 
    %calculate --- L(1,g) 
    L_U_1_4_ = E*L_W_0_4_ - P_U_3_.*E - P_U_3_.*E;  	 
    L_W_1_4_ = Et*L_U_0_4_ - P_W_3_.*Et - P_W_3_.*Et;  
 
    %calculate --- L(g-1,2) 
    L_U_3_2_ = P_U_3_*L_W_0_2_m1_ - P_U_3_.*E; 
    L_W_3_2_ = P_W_3_*L_U_0_2_m1_ - P_W_3_.*Et;  
 
    %calculate --- P(g+1) 
    P_U_5_ = P_U_4_*E - L_U_1_4_ - L_U_3_2_;   
    P_W_5_ = P_U_5_';  
 
    %calculate --- L(0,g+2) 
    L_U_0_6_ = P_U_5_*Et.*Iu; 
    L_W_0_6_ = P_W_5_*E.*Iw; 
 
    if(g_ == 4) 
        Ng2_ = (trace(L_U_0_6_) / 6); 
        %connection of cycles 4 to every variable nodes
        Ng2_per_u_=diag(L_U_0_6_/6);
        %connection of cycles  to every variable nodes
        Ng2_per_p_=diag(L_W_0_6_/ 6);

   
    elseif ((trace(L_U_0_6_) > 0)) 
        g_ = 6; 
        Ng_ = (trace(L_U_0_6_) / 6); 
        %connection of cycles 6 to every variable nodes
        Ng_per_u_=diag(L_U_0_6_/ 6);
        Ng_per_p_=diag(L_W_0_6_/ 6);
 
    end 
 
%Count 8-Cycles------------------------------------------------------------ 
 
    %calculate --- L(2,g) 
    L_U_2_4_ = E*L_W_1_4_ - diag(diag(E*L_W_1_4_)); 
    L_W_2_4_ = Et*L_U_1_4_ - diag(diag(Et*L_U_1_4_)); 
    if(g_ == 4) 
        L_U_2_4_ = L_U_2_4_ - (((P_U_2_).*((P_U_2_)-1).*((P_U_2_)-2))./6 ).*6; 
        L_W_2_4_ = L_W_2_4_ - (((P_W_2_).*((P_W_2_)-1).*((P_W_2_)-2))./6 ).*6; 
    end 
 
    %calculate --- L(g,2) 
    L_U_4_2_ = E*L_W_3_2_ - diag(diag(E*L_W_3_2_)) - (L_U_0_2_m1_*L_U_2_2_); 
    L_W_4_2_ = Et*L_U_3_2_ - diag(diag(Et*L_U_3_2_)) - (L_W_0_2_m1_*L_W_2_2_); 
    if(g_ == 4) 
        L_U_4_2_ = L_U_4_2_ + P_U_2_c2_ + P_U_2_c2_; 
        L_W_4_2_ = L_W_4_2_ + P_W_2_c2_ + P_W_2_c2_;  
    end 
 
    %calculate --- P(g+2) 
    P_U_6_ = P_U_5_*Et - L_U_0_6_ - L_U_2_4_ - L_U_4_2_; 
    P_W_6_ = P_W_5_*E - L_W_0_6_ - L_W_2_4_ - L_W_4_2_; 
 
    %calculate --- L(1,g+2) 
    L_U_1_6_ = E*L_W_0_6_ - (P_U_5_.*E).*2; 
    L_W_1_6_ = Et*L_U_0_6_ - (P_W_5_.*Et).*2; 
    if (g_ == 4) 
        L_U_1_6_ = L_U_1_6_ - (((P_U_3_).*((P_U_3_)-1)./2)).*2.*E ... 
                            + (P_U_2_c2_*E - P_U_3_).*E.*2 ... 
                            + (E*P_W_2_c2_ - P_U_3_).*E.*2; 
 
        L_W_1_6_ = L_W_1_6_ - (((P_W_3_).*((P_W_3_)-1)./2)).*2.*Et ... 
                            + (P_W_2_c2_*Et - P_W_3_).*Et.*2 ... 
                            + (Et*P_U_2_c2_ - P_W_3_).*Et.*2;   
    end 
 
    %calculate --- L(3,g) 
    L_U_3_4_ = L_U_1_4_; 
    if (g_ == 4) 
        L_U_3_4_ = E*L_W_2_4_ - (L_U_0_2_m1_*L_U_1_4_) ... 
                              - (((P_U_3_).*((P_U_3_)-1)./2)).*4.*E ... 
                              + (P_U_2_c2_*E - P_U_3_).*E.*4 ... 
                              + (E*P_W_2_c2_ - P_U_3_).*E.*6;    
    end 
 
    %calculate --- L(g+1,2) 
    L_U_5_2_ = E*L_W_4_2_ - (L_U_0_2_m1_*L_U_3_2_) - P_U_5_.*E; 
    L_W_5_2_ = Et*L_U_4_2_ - (L_W_0_2_m1_*L_W_3_2_) - P_W_5_.*Et; 
    if (g_ == 4)  
        L_U_5_2_ = L_U_5_2_ + P_U_3_.*E + P_U_3_.*E ... 
                                        - (L_U_0_4_*L_U_1_2_) ... 
                                        + (L_U_3_2_.*E)... 
                                        + (P_U_3_*L_W_0_2_m2_).*E ... 
                                        + (P_U_3_*L_W_0_2_m2_).*E ... 
                                        + ((P_U_2_c2_*E) - P_U_3_).*E ... 
                                        + ((P_U_2_c2_*E) - P_U_3_).*E; 
 
        L_W_5_2_ = L_W_5_2_ + P_W_3_.*Et + P_W_3_.*Et ... 
                                         - (L_W_0_4_*L_W_1_2_) ... 
                                         + (L_W_3_2_.*Et)... 
                                         + (P_W_3_*L_U_0_2_m2_).*Et ... 
                                         + (P_W_3_*L_U_0_2_m2_).*Et ... 
                                         + ((P_W_2_c2_*Et) - P_W_3_).*Et ... 
                                         + ((P_W_2_c2_*Et) - P_W_3_).*Et; 
    end 
 
    %calculate --- P(g+3) 
    P_U_7_ = P_U_6_*E - L_U_1_6_ - L_U_3_4_ - L_U_5_2_; 
    P_W_7_ = P_U_7_'; 
 
    %calculate --- L(0,g+4) 
    L_U_0_8_ = P_U_7_*Et.*Iu; 
    L_W_0_8_ = P_W_7_*E.*Iw; 
 
    if(g_ == 4)         
        Ng4_ = (trace(L_U_0_8_) / 8); 
        %connection of cycles to every variable nodes
        Ng4_per_u_=diag(L_U_0_8_/8);
        Ng4_per_p_=diag(L_W_0_8_/8);


        return;         
    elseif(g_ == 6)         
        Ng2_ = (trace(L_U_0_8_) / 8);      
        Ng2_per_u_=diag(L_U_0_8_/8);
        Ng2_per_p_=diag(L_W_0_8_/8);


    elseif ((trace(L_U_0_8_) > 0) && (g_~= 6)) 
        g_ = 8; 
        Ng_ = (trace(L_U_0_8_) / 8);       
        Ng_per_u_=diag(L_U_0_8_/8);
        Ng_per_p_=diag(L_W_0_8_/8);

    end 
     
%Count 10-Cycles for g6 and 12-Cycles for g8------------------------------- 
 
if (g_ > 4) 
     
    % count 10-Cycles for g6----------------------------------------------- 
     
    L_U_2_g_    = E*L_W_1_6_    - diag(diag(E*L_W_1_6_)); 
	L_W_2_g_    = Et*L_U_1_6_   - diag(diag(Et*L_U_1_6_));	 
 
    L_U_g_2_    = E*L_W_5_2_  - diag(diag(E*L_W_5_2_)) ... 
                                + P_U_4_.*P_U_2_ ... 
                                - L_U_0_2_m1_*L_U_4_2_; 
                             
    L_W_g_2_    = Et*L_U_5_2_ - diag(diag(Et*L_U_5_2_))... 
                                + P_W_4_.*P_W_2_ ... 
                                - L_W_0_2_m1_*L_W_4_2_; 
 
    P_U_g2_     = P_U_7_*Et    - L_U_0_8_ - L_U_2_g_ - L_U_g_2_; 
    
    L_U_1_g2_   = E*L_W_0_8_   - P_U_7_.*2.*E; 
     
    L_U_3_g_    = E*L_W_2_g_    - L_U_0_2_m1_*L_U_1_6_ ... 
                                - ((P_U_3_.*(P_U_3_-1).*(P_U_3_-2))).*6;  
     
    L_U_g1_2_	= E*L_W_g_2_    - P_U_7_.*E ... 
                                + L_U_5_2_.*E ... 
                                - L_U_0_2_m1_*L_U_5_2_ ... 
                                - L_U_0_6_*L_U_1_2_ ... 
                                + 2.*E.*P_U_5_ ... 
                                + (P_U_5_*L_W_0_2_m2_).*2.*E; 
                             
    P_U_g3_     = P_U_g2_*E - L_U_1_g2_ - L_U_3_g_ - L_U_g1_2_; 
     
    L_U_0_g4_	= P_U_g3_*Et.*Iu; 
    

    P_W_g2_     = P_W_7_*E - L_W_0_8_ - L_W_2_g_ - L_W_g_2_; 
         
    L_W_1_g2_   = Et*L_U_0_8_ - P_W_7_.*2.*Et; 
         
    L_W_g1_2_   = Et*L_U_g_2_   - P_W_7_.*Et ... 
                                    + L_W_5_2_.*Et ... 
                                    - L_W_0_2_m1_*L_W_5_2_ ... 
                                    - L_W_0_6_*L_W_1_2_ ... 
                                    + P_W_5_.*2.*Et ... 
                                    + (P_W_5_*L_U_0_2_m2_).*2.*Et; 
                                 
    L_W_3_g_    = Et*L_U_2_g_   - L_W_0_2_m1_*L_W_1_6_ ... 
                                    - ((P_W_3_.*(P_W_3_-1).*(P_W_3_-2))./6).*6; 
         
    P_W_g3_     = P_W_g2_*Et - L_W_1_g2_ - L_W_3_g_ - L_W_g1_2_;     
         
    L_W_0_g4_   = P_W_g3_*E.*Iw; 
    
    if (g_ == 6) 
        Ng4_    = trace(L_U_0_g4_)/ (g_+ 4); 
        
        Ng4_per_u_=diag(L_U_0_g4_)/ (g_+ 4);
        Ng4_per_p_=diag(L_W_0_g4_)/ (g_+ 4);

    
        return;   
    else 
 
        % count 12-Cycles for g8 
         
        P_W_g2_     = P_W_7_*E - L_W_0_8_ - L_W_2_g_ - L_W_g_2_; 
         
        L_W_1_g2_   = Et*L_U_0_8_ - P_W_7_.*2.*Et; 
         
        L_W_g1_2_   = Et*L_U_g_2_   - P_W_7_.*Et ... 
                                    + L_W_5_2_.*Et ... 
                                    - L_W_0_2_m1_*L_W_5_2_ ... 
                                    - L_W_0_6_*L_W_1_2_ ... 
                                    + P_W_5_.*2.*Et ... 
                                    + (P_W_5_*L_U_0_2_m2_).*2.*Et; 
                                 
        L_W_3_g_    = Et*L_U_2_g_   - L_W_0_2_m1_*L_W_1_6_ ... 
                                    - ((P_W_3_.*(P_W_3_-1).*(P_W_3_-2))./6).*6; 
         
        P_W_g3_     = P_W_g2_*Et - L_W_1_g2_ - L_W_3_g_ - L_W_g1_2_;     
         
        L_W_0_g4_   = P_W_g3_*E.*Iw; 
         
        P_U_g1_     =  P_U_g3_; 
        %P_W_g1_     = P_W_g3_; 
        L_U_0_g2_   = L_U_0_g4_; 
        L_W_0_g2_   = L_W_0_g4_; 
        L_U_gm2_2_  = L_U_g_2_; 
        L_W_gm2_2_  = L_W_g_2_;  
        L_U_gm1_2_  = L_U_g1_2_; 
        L_W_gm1_2_  = L_W_g1_2_; 
        L_U_1_g_    = L_U_1_g2_; 
        L_W_1_g_    = L_W_1_g2_;   
         
         
        L_U_2_g_    = E*L_W_1_g_    - diag(diag(E*L_W_1_g_)); 
        L_W_2_g_    = Et*L_U_1_g_   - diag(diag(Et*L_U_1_g_));	 
 
        L_U_g_2_    = E*L_W_gm1_2_  - diag(diag(E*L_W_gm1_2_)) ... 
                                    - L_U_0_2_m1_*L_U_gm2_2_ ... 
                                   + P_U_6_.*P_U_2_;                       
                                 
        L_W_g_2_    = Et*L_U_gm1_2_ - diag(diag(Et*L_U_gm1_2_)) ... 
                                    - L_W_0_2_m1_*L_W_gm2_2_ ... 
                                    + P_W_6_.*P_W_2_; 
 
                                 
        P_U_g2_      = P_U_g1_*Et	- L_U_0_g2_ - L_U_2_g_ - L_U_g_2_; 
 
 
        L_U_1_g2_   = E*L_W_0_g2_	- 2.*P_U_g1_.*E;      
         
        L_U_3_g_    = E*L_W_2_g_	- L_U_0_2_m1_*L_U_1_g_;                             
                                     
        L_U_g1_2_   = E*L_W_g_2_	- P_U_g1_.*E ... 
                                    + L_U_gm1_2_.*E ... 
                                    - L_U_0_2_m1_*L_U_gm1_2_ ... 
                                    - L_U_0_8_*L_U_1_2_ ... 
                                    + 2.*P_U_7_.*E ... 
                                    + (P_U_7_*L_W_0_2_m2_).*2.*E;                        
                                     
                                                              
        P_U_g3_     = P_U_g2_*E - L_U_1_g2_ - L_U_3_g_ - L_U_g1_2_; 
 
        L_U_0_g4_	= P_U_g3_*Et.*Iu; 
 
        Ng2_        = trace(L_U_0_g2_) / (g_+ 2); 
        Ng2_per_u_=diag(L_U_0_g2_/(g_+ 2));
        Ng2_per_p_=diag(L_W_0_g2_/(g_+ 2));
        
         
        Ng4_        = trace(L_U_0_g4_) / (g_+ 4); 
        Ng4_per_u_=diag(L_U_0_g4_ / (g_+ 4));
        Ng4_per_p_=diag(L_W_0_g4_ / (g_+ 4));
        %statistics of short cycles
        %Ng_mean_VN = mean(Ng_per_p_);
        %Ng2_mean_VN = mean(Ng2_per_p_);
        %Ng4_mean_VN = mean(Ng4_per_p_);
        %Ng_mean_CN = mean(Ng_per_u_);
        %Ng2_mean_CN = mean(Ng2_per_u_);
        %Ng4_mean_CN = mean(Ng4_per_u_);
        %Ng_std_VN_= std(Ng_per_p_);
        %Ng2_std_VN = std(Ng2_per_p_);
        %Ng4_std_VN = std(Ng4_per_p_);
        %Ng_std_CN=std(Ng_per_u_);
        %Ng2_std_CN=std(Ng2_per_u_);
        %Ng4_std_CN=std(Ng4_per_u_);
    
    
    
    
        return; 

    end 


end 