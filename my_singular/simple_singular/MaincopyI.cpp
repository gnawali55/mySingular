#include <singular/Singular/libsingular.h>
#include <iostream>
#include <cstring>  // For strdup



lists pOLY_List(poly f0) 
{
    if (f0 == NULL) return NULL; // Handle NULL input

    int r = pLength(f0);
    lists S = (lists)omAlloc0Bin(slists_bin);
    S->Init(r);

    for (int k = 0; k < r; k++) 
    {
        S->m[k].rtyp = POLY_CMD;  // Ensure correct type

        poly headTerm = pHead(f0);
        if (headTerm != NULL) 
        {
            S->m[k].data = (void*) p_Copy(headTerm, currRing); // Deep copy
        } 
        else 
        {
            std::cout << "Warning: pHead(f0) returned NULL at index " << k << std::endl;
        }

        f0 = pNext(f0);
        if (f0 == NULL) break;
    }

    return S;
}


matrix lcm_mod(ideal G) { //ideal G is Singular module

    int a=0;
    int b=0;
    int i=0;
    int j=0;
    ideal G_copy=idCopy(G);
    int r = IDELEMS(G_copy);
    matrix l=mpNew(r,r);
    poly s10=NULL;
    poly t10=NULL;
    for (a = 0; a < r; a++) {
        for (b = 0; b < r; b++) {
             //std::cout << "G->m[a]>: " << pString(G->m[a]) << std::endl;
             //std::cout << "G->m[b]>: " << pString(G->m[b]) << std::endl;
           
            //i = leadexp(G[a])[nvars(basering) + 1];
            i= p_GetComp(G->m[a],currRing);
            //j = leadexp(G[b])[nvars(basering) + 1];
            j= p_GetComp(G->m[b],currRing);
            
   
            s10 = pHead(G->m[a]); //Should be leadmonomial
            pSetComp(s10,0);
            pSetmComp(s10);
            t10 = pHead(G->m[b]);//Should be leadmonomial
            pSetComp(t10,0);
            pSetmComp(t10);
            poly lcm_poly =  p_Lcm(s10, t10, currRing);
            pSetCoeff0(lcm_poly,nInit(1));

            if (i == j) {
               // l[a, b] = lcm(leadmonomial(G[a]), leadmonomial(G[b])) / lead(t10);
                
                MATELEM(l, a, b) = pp_Divide(lcm_poly, t10, currRing);
                //  std::cout << "m[a.b] in lcm_mod: " << pString(MATELEM(l, a, b)) << std::endl;

            } else {
                // If i is not equal to j, set l[a, b] to 0
                MATELEM(l,a,b)= NULL;
            }
        }
    }
    

    return l;
}





ideal leadSyz(ideal f) {
    int a = 0, b = 0, j = 0, k = 0;  // Loop variables
    poly s = NULL;  // Temporary polynomial to store elements of L
    poly t = NULL;  // Temporary polynomial to store the current syzygy candidate
    ideal f_copy = idCopy(f);  // Copy the input ideal to avoid modifying it
    int r = IDELEMS(f_copy);  // Number of elements in the ideal
    ideal L = idInit(0, 1);  // Initialize an empty ideal (syzygy module)
    matrix M = mpNew(r, r);  // Matrix to store LCM-based computations

    // Fill the matrix M with LCM computations
    for (a = 0; a < r; a++) {
        for (b = 0; b < r; b++) {
            // Compute LCM of leading monomials of f[a] and f[b]
            poly lcm = p_Lcm(pHead(f_copy->m[b]), pHead(f_copy->m[a]), currRing);
            pSetCoeff0(lcm, nInit(1));  // Normalize the LCM (set coefficient to 1)
            MATELEM(M, a, b) = pp_Divide(lcm, pHead(f_copy->m[b]), currRing); // Store the quotient
        }
    }

    int cc = 0;  // Counter for the number of generators in L

    // Iterate through pairs of indices to construct syzygies
    for (int i = 1; i < r; i++) {
        for (j = 0; j < i; j++) {
            // Generate the initial syzygy candidate from matrix M
            poly t0 = pCopy(MATELEM(M, j, i));
            p_SetComp(t0, i + 1, currRing);  // Assign the component index
            p_SetmComp(t0, currRing);       // Normalize the component
            t = pCopy(t0);  // Copy t0 to t as the current syzygy candidate

            // Check divisibility conditions for t against the elements in L
            for (k = 0; k < IDELEMS(L); k++) {
                s = (poly)L->m[k];  // Retrieve the k-th generator of L
                bool c = p_DivisibleBy(pHead(s), t, currRing);  // Check if s divides t
                if (c == TRUE) {
                    // If s divides t, discard t
                    t = NULL;
                    break;
                } else {
                    // Check if t divides s
                    bool d = p_DivisibleBy(pHead(t), s, currRing);
                    if (d == TRUE) {
                        // If t divides s, remove s from L
                        // std::cout << "t divides L[" << k << "], removing L[" << k << "]." << std::endl;
                        ideal tmp = id_Delete_Pos(L, k, currRing);  // Remove s from L
                        L = NULL;  // Set L to NULL before reassigning
                        idDelete(&L);  // Free the old L
                        L = tmp;  // Assign the updated L
                        k--;  // Adjust index to account for the removed element
                        cc--;  // Decrease the counter
                    }
                }
            }

            // If t survives all checks, add it to L
            if (t != NULL) {
                if (cc >= IDELEMS(L)) {
                    // Resize L if necessary
                    ideal tmpL = idInit(cc + 1, 1);
                    for (int i = 0; i < cc; i++) {
                        tmpL->m[i] = pCopy(L->m[i]);  // Copy elements from L to tmpL
                        L->m[i] = NULL;  // Clear the old L
                    }
                    idDelete(&L);  // Free the old L
                    L = tmpL; 
                     // Update L to point to the resized ideal
                }
                L->m[cc] = pCopy(t);  // Add t to L
                cc++;  // Increment the counter
                p_Delete(&t, currRing);  // Free t after adding it to L
            }
        }
    }

    // Debug output: Print the final size and contents of L
    std::cout << "Final first leadsyz size: " << IDELEMS(L) << std::endl;
    for (int k = 0; k < cc; k++) {
        std::cout << "Generator " << k << ": " << pString((poly)L->m[k]) << std::endl;
    }

 

    return L;  // Return the computed syzygy module
}






ideal Sec_leadSyz(ideal f0) {
    int r = IDELEMS(f0);  // Get the number of elements in the ideal f0
    poly s = NULL;  // Polynomial s is a singular vector
    poly t = NULL;  // Polynomial t is a singular vector
    int cc = 0; // Counter for elements in L
    
    // Initialize ideal L with initial size 0
    ideal L = idInit(0, 1);
    ideal f_copy=idCopy(f0);
    // Create a matrix M using lcm_mod for the input ideal f0
    matrix M = lcm_mod(f_copy);  // Ensure lcm_mod returns a valid matrix

    // Loop through pairs (i, j) in the matrix
    for (int i = 1; i < r; i++) {
        for (int j = 0; j < i; j++) {
            // Fetch the matrix element at (j, i)
            poly t0 = pCopy(MATELEM(M, j, i));

            if (t0 != NULL) {
                // Set the component and multigrade component for t0
                p_SetComp(t0, i + 1, currRing);
                p_SetmComp(t0, currRing);

                // Copy t0 into t
                t = pCopy(t0);  
            }

            // Ensure L is not NULL before accessing it
            for (int k = 0; k < IDELEMS(L); k++) { 
                // Fetch the k-th element of L (s = L[k])
                s = (poly)L->m[k];  

                // Ensure both s and t are not NULL before checking divisibility
                if (s != NULL && t != NULL) {
                    // Check if s divides t
                    if (p_DivisibleBy(pHead(s), t, currRing)) {
                        // If s divides t, set t to NULL and break out of the loop
                        t = NULL; 
                        break;
                    } 
                    // Check if t divides s
                    else if (p_DivisibleBy(t, pHead(s), currRing)) {
                        // Log the removal for debugging purposes
                        // std::cout << "Removing s =: " << pString(s) << ": t=" << pString(t) << std::endl;

                        // Remove s from L using id_Delete_Pos
                        ideal tmp = id_Delete_Pos(L, k, currRing);
                        
                        // Set L to NULL, delete the old L, and assign the new ideal
                        idDelete(&L);  // Delete old L
                        L = tmp;  // Assign new ideal

                        // Adjust indexing and counter after deletion
                        k--;  // Reindex to avoid skipping elements
                        cc--; // Decrement the counter after the deletion
                    }
                }
            }

            // If t is not NULL, add it to L
            if (t != NULL) {
                // Resize L if necessary
                if (cc >= IDELEMS(L)) {
                    // Create a temporary ideal tmpL with space for one more element
                    ideal tmpL = idInit(cc + 1, 1);

                    // Copy elements from L to tmpL
                    for (int i = 0; i < cc; i++) {
                        tmpL->m[i] = pCopy(L->m[i]);  // Copy elements
                        L->m[i] = NULL;  // Clear L after copying
                    }

                    // Delete old L and assign the resized tmpL to L
                    idDelete(&L);  // Delete old L
                    L = tmpL;  // Assign resized tmpL to L
                }

                // Add t to L at the next available position
                L->m[cc] = t;
                cc++;  // Increment counter after adding t
            } 
        }
    }

      //  Debug output
    std::cout << "Final second-Leadsyz  size: " << IDELEMS(L) << std::endl;
    // for (int k = 0; k < cc; k++) {
    //     std::cout << "Generator " << k << ": " << pString((poly)L->m[k]) << std::endl;
    // }

    // // return (L);
    return L;
}



lists aLL_LEAD(ideal f) {
    // Allocate memory for the lists structure
    lists J = (lists)omAlloc0Bin(slists_bin);
    J->Init(2); // Initialize the list with two elements
    ideal f_copy = f;
    
    // Initialize the first two elements
    J->m[0].rtyp = IDEAL_CMD;  
    J->m[0].data = f_copy;          
    
    int n = rVar(currRing); // Get the number of variables in the current ring
    ideal F = leadSyz(f_copy);
    int g=IDELEMS(F);
   
    for(int k=0; k<g;k++){
      std::cout << "First_LeadSyz :at" <<""<<k << ": " << pString((poly)F->m[k]) << std::endl;
    }
    ideal F_copy =F;
   
    J->m[1].rtyp = MODUL_CMD;  
    J->m[1].data = F_copy;

    int cc = 2; // Current count of elements in the list

    // Debugging output to print the number of variables
    // std::cout << "Number of variables: " << n << std::endl;

    // Iterate through the variables in the ring to compute subsequent syzygy ideals
    for (int t = 0; t < n; t++) {
        // std::cout << "Processing variable: " << t + 1 << std::endl;

        
        ideal temp = Sec_leadSyz(F_copy);
        bool b=idIs0(temp);
       if(b==FALSE){
        // idDelete(&F_copy);
         F_copy = temp; 
         temp=NULL;
       }
    else
    {
            std::cout << "In break: the syzygy ideal is zero or empty." << std::endl;
            break;  // Exit the loop
    }

        // Resize the list if necessary to accommodate new elements
        if (cc >= lSize(J)+1) {
            int newSize = cc+1; // Increase size by 1
            lists tmpL = (lists)omAlloc0Bin(slists_bin); // Allocate a new list
            tmpL->Init(newSize); // Initialize the new list with the updated size

            // Copy elements from the old list `J` to the new list `tmpL`
            for (int i = 0; i < cc; i++) {
                //  std::cout << "variable:i " << i<< std::endl;
                tmpL->m[i].rtyp = J->m[i].rtyp;
                tmpL->m[i].data = J->m[i].data;
            }

            // Free memory associated with the old list
            omFreeBin(J, slists_bin);

            // Assign the resized list to `J`
            J = tmpL;
        }
     std::cout << "Counter in the loop:= " << cc << std::endl;
        // Append the new syzygy ideal to the list
        J->m[cc].rtyp = MODUL_CMD;
        J->m[cc].data = F_copy;
            cc++;
        
    }
       for(int k=0;k<cc;k++){
      ideal l=(ideal)J->m[k].Data();
      for(int s=0; s <IDELEMS(l);s++){
      std::cout << "Sch FRame as list " << k << ": " << pString((poly)l->m[s]) << std::endl;
      }
 
    }
    return J;
}int main() {
    siInit((char *)"/home/santosh/singular-gpispace/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/singular-4.4.0p2-syrkttc4im2j3tzob5jykruuxnushksj/lib/libSingular.so");

    // Define the variables
    char **n = (char**)omalloc(4 * sizeof(char*));
    n[0] = omStrDup("w");
    n[1] = omStrDup("x");
    n[2] = omStrDup("y");
    n[3] = omStrDup("z");

    // Define the ring (dp,c) ordering
    ring R = rDefault(0, 4, n);  // 0 means coefficient field is ℚ
    rChangeCurrRing(R);

    // Define the polynomials for the ideal
    poly f1 = p_ISet(1, R);
    pSetExp(f1, 1, 2); // w^2

    pSetm(f1);
    

   // Create the polynomial -1 x*z
  poly p2 = p_ISet(-1, R);
  pSetExp(p2, 2, 1);
  pSetExp(p2, 4, 1);
 
  pSetm(p2);

   f1 = p_Add_q(f1, p2, R);
  p2 = NULL;
  pWrite(f1); //w^2-x*z



// Define f2 = w*x - y*z
poly f2 = p_ISet(1, R);
pSetExp(f2, 1, 1); // w
pSetExp(f2, 2, 1); // x
pSetm(f2);

poly p3 = p_ISet(-1, R); // -y*z
pSetExp(p3, 3, 1); // y
pSetExp(p3, 4, 1); // z
pSetm(p3);

f2 = p_Add_q(f2, p3, R);
 p3 = NULL;
pWrite(f2);
printf("\n");

// Define f3 = x^2 - w*y
poly f3 = p_ISet(1, R);
pSetExp(f3, 2, 2); // x^2
pSetm(f3);

poly p4 = p_ISet(-1, R); // -w*y
pSetExp(p4, 1, 1); // w
pSetExp(p4, 3, 1); // y
pSetm(p4);

f3 = p_Add_q(f3, p4, R);
 p4 = NULL;
pWrite(f3);
printf("\n");

// Define f4 = x*y - z^2
poly f4 = p_ISet(1, R);
pSetExp(f4, 2, 1); // x
pSetExp(f4, 3, 1); // y
pSetm(f4);

poly p5 = p_ISet(-1, R); // -z^2
pSetExp(p5, 4, 2); // z^2
pSetm(p5);

f4 = p_Add_q(f4, p5, R);
 p5 = NULL;
pWrite(f4);
printf("\n");

// Define f5 = y^2 - w*z
poly f5 = p_ISet(1, R);
pSetExp(f5, 3, 2); // y^2
pSetm(f5);

poly p6 = p_ISet(-1, R); // -w*z
pSetExp(p6, 1, 1); // w
pSetExp(p6, 4, 1); // z
pSetm(p6);

f5 = p_Add_q(f5, p6, R);
 p6 = NULL;
pWrite(f5);
printf("\n");

// Construct the ideal J
ideal J = idInit(5, 1);
J->m[0] = f1;
J->m[1] = f2;
J->m[2] = f3;
J->m[3] = f4;
J->m[4] = f5;



    // Run the aLL_LEAD function
    lists syzygyList = aLL_LEAD(J);

    // Print the results
    std::cout << "Computed Syzygies: " << std::endl;
    for (int i = 0; i < lSize(syzygyList); i++) {
        ideal result = (ideal)syzygyList->m[i].data;
        std::cout << "Level " << i << " Syzygies: " << std::endl;
        for (int j = 0; j < IDELEMS(result); j++) {
            std::cout << "  " << pString((poly)result->m[j]) << std::endl;
        }
    }

    // Cleanup
    idDelete(&J);
    rKill(R);
    return 0;
}
