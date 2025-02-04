#include <singular/Singular/libsingular.h>
#include <iostream>
#include <cstring>  // For strdup

matrix lcm_mod(ideal G) {
    int a = 0, b = 0, i = 0, j = 0;
    ideal G_copy = idCopy(G);
    int r = IDELEMS(G_copy);
    for (int k = 0; k < r; k++) {
        std::cout << "ideal G: " << pString((poly)G->m[k]) << std::endl;
    }
    std::cout << "r in lcm_mod:= " << r << std::endl;
    omUpdateInfo();
    std::cout << "used mem: " << om_Info.UsedBytes << std::endl;

    matrix l = mpNew(r, r);
    std::cout << "pointer l: " << (long)l << std::endl;
    std::cout << "Row of l: " << MATROWS(l) << std::endl;
    std::cout << "Cols of l: " << MATCOLS(l) << std::endl;

    poly s10 = NULL;
    poly t10 = NULL;

    for (a = 1; a < r+1; a++) {
        for (b = 1; b < r+1; b++) {
            i = p_GetComp(G->m[a-1], currRing);
            j = p_GetComp(G->m[b-1], currRing);

            s10 = pHead(G->m[a-1]);
            pSetComp(s10, 0);
            pSetmComp(s10);

            t10 = pHead(G->m[b-1]);
            pSetComp(t10, 0);
            pSetmComp(t10);

            poly lcm_poly = p_Lcm(s10, t10, currRing);
            pSetCoeff0(lcm_poly, nInit(1));

            if (i == j) {
                MATELEM(l, a, b) = pp_Divide(lcm_poly, t10, currRing);
            } else {
                MATELEM(l, a, b) = NULL;
            }
        }
    }

    p_Delete(&s10, currRing);
    p_Delete(&t10, currRing);

    return l;
}

ideal leadSyz(ideal f) {
    int a = 0, b = 0, j = 0, k = 0;
    poly s = NULL, t = NULL;
    ideal f_copy = idCopy(f);
    int r = IDELEMS(f_copy);
    ideal L = idInit(0, 1);
    matrix M = mpNew(r, r);

    for (a = 1; a < r+1; a++) {
        for (b = 1; b < r+1; b++) {
            poly lcm = p_Lcm(pHead(f_copy->m[b-1]), pHead(f_copy->m[a-1]), currRing);
            pSetCoeff0(lcm, nInit(1));
            MATELEM(M, a, b) = pp_Divide(lcm, pHead(f_copy->m[b-1]), currRing);
        }
    }

    int cc = 0;
    for (int i = 2; i < r+1; i++) {
        for (j = 1; j < i; j++) {
            poly t0 = pCopy(MATELEM(M, j, i));
            p_SetComp(t0, i, currRing);
            p_SetmComp(t0, currRing);
            t = pCopy(t0);

            for (k = 0; k < IDELEMS(L); k++) {
                s = (poly)L->m[k];
                bool c = p_DivisibleBy(pHead(s), t, currRing);
                if (c == TRUE) {
                    t = NULL;
                    break;
                } else {
                    bool d = p_DivisibleBy(pHead(t), s, currRing);
                    if (d == TRUE) {
                        ideal tmp = id_Delete_Pos(L, k, currRing);
                        idDelete(&L);
                        L = tmp;
                        k--;
                        cc--;
                    }
                }
            }

            if (t != NULL) {
                if (cc >= IDELEMS(L)) {
                    ideal tmpL = idInit(cc + 1, 1);
                    for (int i = 0; i < cc; i++) {
                        tmpL->m[i] = pCopy(L->m[i]);
                        L->m[i] = NULL;
                    }
                    idDelete(&L);
                    L = tmpL;
                }
                L->m[cc] = pCopy(t);
                cc++;
                p_Delete(&t, currRing);
            }
        }
    }

    std::cout << "Final first leadsyz size: " << IDELEMS(L) << std::endl;

    return L;
}

ideal Sec_leadSyz(ideal f0) {
    int r = IDELEMS(f0);
    poly s = NULL, t = NULL;
    int cc = 0;
    ideal L = idInit(0, 1);
    ideal f_copy = idCopy(f0);
    matrix M = lcm_mod(f_copy);

    for (int i = 2; i < r+1; i++) {
        for (int j = 1; j < i; j++) {
            poly t0 = MATELEM(M, j, i);
            if (t0 != NULL) {
                p_SetComp(t0, i, currRing);
                p_SetmComp(t0, currRing);
                t = pCopy(t0);
            }

            for (int k = 0; k < IDELEMS(L); k++) {
                s = (poly)L->m[k];
                if (s != NULL && t != NULL) {
                    if (p_DivisibleBy(pHead(s), t, currRing)) {
                        t = NULL;
                        break;
                    } else if (p_DivisibleBy(t, pHead(s), currRing)) {
                        ideal tmp = id_Delete_Pos(L, k, currRing);
                        idDelete(&L);
                        L = tmp;
                        k--;
                        cc--;
                    }
                }
            }

            if (t != NULL) {
                if (cc >= IDELEMS(L)) {
                    ideal tmpL = idInit(cc + 1, 1);
                    for (int i = 0; i < cc; i++) {
                        tmpL->m[i] = pCopy(L->m[i]);
                        L->m[i] = NULL;
                    }
                    idDelete(&L);
                    L = tmpL;
                }
                L->m[cc] = t;
                cc++;
            }
        }
    }

    idDelete((ideal*)&M);

    std::cout << "Final second-Leadsyz size: " << IDELEMS(L) << std::endl;

    return L;
}

lists aLL_LEAD(ideal f) {
    lists J = (lists)omAlloc0Bin(slists_bin);
    J->Init(2);
    ideal f_copy = idCopy(f);
    int n = rVar(currRing);
    ideal F = leadSyz(f_copy);
    int g = IDELEMS(F);
    idDelete(&f);
    ideal F_copy = idCopy(F);

    J->m[0].rtyp = IDEAL_CMD;
    J->m[0].data = f_copy;

    J->m[1].rtyp = MODUL_CMD;
    J->m[1].data = F_copy;

    int cc = 2;

    for (int t = 0; t < n; t++) {
        ideal temp = Sec_leadSyz(F_copy);

        bool b = idIs0(temp);
        if (b == FALSE) {
            F_copy = idCopy(temp);
            temp = NULL;
        } else {
            break;
        }

        if (cc >= lSize(J) + 1) {
            int newSize = cc + 1;
            lists tmpL = (lists)omAlloc0Bin(slists_bin);
            tmpL->Init(newSize);

            for (int i = 0; i < cc; i++) {
                tmpL->m[i].rtyp = J->m[i].rtyp;
                tmpL->m[i].data = J->m[i].data;
            }

            omFreeBin(J, slists_bin);
            J = tmpL;
        }

        J->m[cc].rtyp = MODUL_CMD;
        J->m[cc].data = F_copy;
        cc++;
    }

    return J;
}



int main() {
    siInit((char *)"/home/santosh/singular-gpispace/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/singular-4.4.0p2-syrkttc4im2j3tzob5jykruuxnushksj/lib/libSingular.so");

    // // Define the variables
    // char **n = (char**)omalloc(4 * sizeof(char*));
    // n[0] = omStrDup("w");
    // n[1] = omStrDup("x");
    // n[2] = omStrDup("y");
    // n[3] = omStrDup("z");
    // // n[] = omStrDup("z");
    // rRingOrder_t* order=(rRingOrder_t*)omAlloc0(3*sizeof(rRingOrder_t));
    // order[0]=ringorder_dp;
    // order[1]=ringorder_c;
    // int* block0=(int*)omAlloc(3*sizeof(int));
    // block0[0]=1;
    // int* block1=(int*)omAlloc0(3*sizeof(int));
    // block1[0]=4;

    // // Define the ring (dp,c) ordering
    // ring R = rDefault(0, 4, n,3,order,block0,block1);  // 0 means coefficient field is ℚ
    // rChangeCurrRing(R);


char **n = (char**)omalloc(5 * sizeof(char*)); // Allocate space for 5 variables
n[0] = omStrDup("w");
n[1] = omStrDup("x");
n[2] = omStrDup("y");
n[3] = omStrDup("z");
n[4] = omStrDup("u"); // Add the new variable "u"

 rRingOrder_t* order=(rRingOrder_t*)omAlloc0(3*sizeof(rRingOrder_t));
    order[0]=ringorder_dp;
    order[1]=ringorder_c;
    int* block0=(int*)omAlloc(3*sizeof(int));
    block0[0]=1;
    int* block1=(int*)omAlloc0(3*sizeof(int));
    block1[0]=5;

// Define the ring (dp,c) ordering with 5 variables
ring R = rDefault(0, 5, n, 3, order, block0, block1);  // Now using 5 instead of 4


rChangeCurrRing(R);




    // Define the polynomials for the ideal
    poly f1 = p_ISet(1, R);
    pSetExp(f1, 1, 1); // w

    pSetm(f1);
    

 
  pWrite(f1); 
   



// Define f2 = w*x - y*z
poly f2 = p_ISet(1, R);
pSetExp(f2, 2, 1); // x

pSetm(f2);


pWrite(f2);
printf("\n");

// Define f3 = x^2 - w*y
poly f3 = p_ISet(1, R);
pSetExp(f3, 3, 1); // y
pSetm(f3);

pWrite(f3);
printf("\n");


poly f4 = p_ISet(1, R);
pSetExp(f4, 4, 1); // z

pSetm(f4);

pWrite(f4);
printf("\n");

poly f5 = p_ISet(1, R);
pSetExp(f5, 5, 1); // z

pSetm(f5);

pWrite(f5);
printf("\n");

// Construct the ideal J
ideal J = idInit(5, 1);
J->m[0] = f1;
J->m[1] = f2;
J->m[2] = f3;
J->m[3] = f4;
J->m[4]=f5;



    // Run the aLL_LEAD function
    lists syzygyList = aLL_LEAD(J);

    // Print the results
    std::cout << "Computed Syzygies: " << std::endl;
    for (int i = 0; i < lSize(syzygyList)+1; i++) {
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

