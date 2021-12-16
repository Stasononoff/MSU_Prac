def U_operate(U, coefs, axis):
    if ((len(U) != 2**len(axis)) | (len(coefs) < len(U)) | (max(axis) > len(U))):
        print('Некорректно задано преобразование')
        return 0
    a = coefs
    b = [0]*len(a)
    
    ax = axis[0]
    
    N = len(U)  # размерность преобразования 
    K = len(coefs) # размерность состояния системы
    
    axis_comb = get_combinations(axis)
    old_ind_set = {K+1}
    
    for i in range(K):
        
        zero_index = i&(K-1)
            
        for ax in axis:
            zero_index = zero_index^(1<<ax)
            
        if zero_index in old_ind_set:
            continue
            
        old_ind_set.add(zero_index)

                
        index = zero_index
        
        
        b[index] += a[zero_index]*U[0][0]
        
#         print(axis_comb)
        u1 = 0
        for ax_list in axis_comb:
            u1 += 1
            m = 0
            for ax in ax_list:
                m += 1<<ax
#             print(index, zero_index, m, zero_index^m)
            b[index] += a[zero_index^m]*U[0][u1]

#         print(index,  a[index], a[index^(m<<ax)])
        u0 = 0
        for ax_list in axis_comb:
            u0 += 1
            r = 0
            for ax in ax_list:
                r += 1<<ax
                
            index = zero_index^r
#             print(index)
            old_ind_set.add(index)
            
            b[index] += a[zero_index]*U[u0][0]
            print('b[{index}]')
            
            u1 = 0
            for ax_list in axis_comb:
                u1 += 1
                m = 0
                for ax in ax_list:
                    m += 1<<ax
                b[index] += a[zero_index^m]*U[u0][u1]
                print()
#         print(old_ind_set)

    #         print(index,  a[index], a[index^(m<<ax)])
        
    return b
        

