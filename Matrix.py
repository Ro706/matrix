import numpy as np  
import math as m 
import time 

print("Want To Try Your Matrix Data Press Y or y\n Else press any key")
g=input(">")
if g == "y"or g == "Y":
    print("start putting Data Line Wise")
    print("")
    user = []
    for i in range(9):
        try:
            row= float(input(f'{i+1}>>> '))
            user.append(row)
        except:
            print('Wrong Input')
            row= float(input(f'{i+1}>>> '))
            user.append(row)
    print('')
    a =np.array([[user[0],user[1],user[2]],[user[3],user[4],user[5]],[user[6],user[7],user[8]]])
else:
    a = np.array([[5,-2,1],[3,1,-4],[6,0,-3]])
    print('')
def print_mat(a):
    print('Given Matrix','\n')
    print('A = ')
    return np.array(a)
def Shape(a):
    print('')
    return 'Given Matrix shape '+str(np.shape(a))
def Dim(a):
    print('')
    return a.ndim
def Determinants_first_value_given(a):
    # First Determinants_first_value_given Take num1[from Matrx]xwith all data 
    # and return the math calculation of sum of the Matrx
    First = a[0][0]*(a[1][1]*a[2][-1]-a[1][-1]*a[-1][1])
    pr_nt1= f'{a[0][0]}{a[1][1]} x {a[2][-1]} - {a[1][-1]} x {a[-1][1]}'
    return  First#,pr_nt1
def F_1(a):
    pr_nt1= f'{a[0][0]}{a[1][1]} x {a[2][-1]} - {a[1][-1]} x {a[-1][1]}'
    pr_nt2 = f'-{a[0][1]}{a[1][0]} x {a[-1][-1]} - {a[1][-1]} x {a[-1][0]}'
    pr_nt3 = f'{a[0][-1]} {a[1][0]}x{a[-1][1]}-{a[1][1]} x {a[-1][0]}'
    return pr_nt1+"-"+pr_nt2+pr_nt3
    
def Determinants_Second_value_given(a):
    # Second Determinants_Second_value_given Take num1[from Matrx]xwith all data 
    # and return the math calculation of sum of the Matrx
    Second = a[0][1]*(a[1][0]*a[-1][-1]-a[1][-1]*a[-1][0])
    pr_nt2 = f'-{a[0][1]}{a[1][0]} x {a[-1][-1]} - {a[1][-1]} x {a[-1][0]}'
    return Second#,pr_nt2
def Determinants_Third_value_given(n):
    # Third Determinants_Thired_value_given Take num1[from Matrx]xwith all data 
    # and return the math calculation of sum of the Matrx
    Third = (a[0][-1]*(a[1][0]*a[-1][1]-a[1][1]*a[-1][0]))
    pr_nt3 = f'{a[0][-1]} {a[1][0]}x{a[-1][1]}-{a[1][1]} x {a[-1][0]}'
    return Third#,pr_nt3
def Matrix_data_convert_ins_data(a):
    print('')
    # Take The data and find The A-1 of the given Data or NUM of the matrix 
    # IM using a inbuild method called a alge of the np data 
    # np.linalg.inv(matrix)
    print('Inverse of A =')
    print('')
    return np.linalg.inv(a)
def s(i):
    return str(i)
def Minors_of_given_matrix_Cofactors(a):
    # The Minors_of_given_matrix is a fun thats return The minors or A[nxn-mxm] of given data 
    # The use of this To get a 1/|A| and comp between A-1 == 0 or not == 0 
    # minor is a data consider or taken by The short Ther Mnxn
    print('')
    t(4)
    print('A1 =','\n',np.array(a))
    m11 = (a[1][1]*a[2][-1]-a[1][-1]*a[-1][1])
    m12 = a[1][0]*a[-1][-1]-a[1][-1]*a[-1][0]
    m13 = (a[1][0]*a[-1][1]-a[1][1]*a[-1][0])
    m21 = a[0][1]*a[-1][-1]-a[0][-1]*a[-1][1]
    m22 = a[0][0]*a[-1][-1]-a[0][-1]*a[-1][0]
    m23 = a[0][0]*a[-1][1]-a[0][1]*a[-1][-1]
    m31 = a[0][1]*a[1][-1]-a[0][-1]*a[1][1]
    m32 = a[0][0]*a[1][-1]-a[0][-1]*a[1][0]
    m33 = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    print('')
        #Now Finding The Cofactors_of_given_matrix
    #Cofactors_of_given_matrix can Find by sq of Minors
    c11 = m.pow((-1),1+1)*m11
    c12 = m.pow((-1),1+2)*m12
    c13 = m.pow((-1),1+3)*m13
    c21 = m.pow((-1),2+1)*m21 
    c22 = m.pow((-1),2+2)*m22 
    c23 = m.pow((-1),2+3)*m23
    c31 = m.pow((-1),3+1)*m31 
    c32 = m.pow((-1),3+2)*m32 
    c33 = m.pow((-1),3+3)*m33
    t(4)
    print('Minors_data = ','\nM11 = '+str(m11),'M12 = '+str(m12),'M13 = '+str(m13),'\nM21 = '+str(m21),'M22 = '+s(m22),'M23 = '+s(m23),'\nM31 = '+s(m31),'M32 = '+s(m32),'M33 = '+s(m33))
    t(3)
    print('')
    print('Cofactors_data =','\nc11 ='+str(c11),'c12 ='+str(c12),'c13 ='+str(c13),'\nc21 ='+str(c21),'c22 ='+str(c22),'c23 ='+str(c23),'\nc31 ='+str(c31),'c32 ='+str(c32),'c33 ='+str(c33))
    t(4)
    print('')
    print('Minors =','\n',np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]))
    print('')
    # print('Cofactors','M11 ='+str(c11),'M12 ='+str(c12),'M13 ='+str(c13),'M21 ='+str(c21),'M22 ='+str(c22),'M23 ='+str(c23),'M31 ='+str(c31),'M32 ='+str(c32),'M33 ='+str(c33))
    t(4)
    print('')
    print('Cofactors =','\n',(np.array([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]],dtype='float')))
    # Adjoint Of the Matrix 
    # TO Adjoint r--.> c row = col 
    print('')
    t(4)
    print('Adjoint = ','\n',np.array([[c11,c21,c31],[c12,c22,c32],[c13,c23,c33]],dtype='float'))
    # find the Invers exixt or not 
    Invers = False
    Invers_find = np.linalg.det(a)
    print('')
    print('|A|=',Invers_find)
    print('')
    if Invers_find >0:
        Invers = True
        print('|A| not 0 So Inverse Existing')
        t(3)
        print(Matrix_data_convert_ins_data(a))
    elif Invers_find <0:
        print('|A| not 0 So Inverse Existing')
        Invers = True
        t(3)
        print(Matrix_data_convert_ins_data(a))
    else:
        print('|A| is = 0 ')
    return
# print(Determinants_first_value_given(a)[0]- Determinants_Second_value_given(a)[0]+Determinants_Third_value_given(a)[0])

# print(Minors_of_given_matrix_Cofactors(a))
# print(Cofactors_of_given_matrix(a))
def t(n):
    return time.sleep(n)
def MAIN_FUN(a):
    print("Sadly, in the next 2 minutes for  Matrix Solver.")
    print('')
    t(4)
    print(print_mat(a))
    t(3)
    print(Shape(a))
    print('The Give Matrix dimensions '+str(Dim(a))+'D')
    t(4)
    print('')
    print("Lets Solve The Matrix .. ... ")
    t(4)
    print(Minors_of_given_matrix_Cofactors(a))
    print('')
    t(3)
    print('We Get ALL Data We need To Solve The Matrix__')
    print(print_mat(a))
    print('')
    t(4)
    print('det(A) =',F_1(a))
    print('')
    t(3)
    print('det(A)=',Determinants_first_value_given(a)- Determinants_Second_value_given(a)+Determinants_Third_value_given(a))
    print('')
    t(3)
    return f'The det(A) value Is {Determinants_first_value_given(a)- Determinants_Second_value_given(a)+Determinants_Third_value_given(a)}'
    
print(MAIN_FUN(a))

print("")
