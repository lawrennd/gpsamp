load results/bootstrap

drosPrintBootstrapResults(res1(:, :, 1:3), 1:4, {'MAP', 'ML', 'Reg', 'Inf'}, T1(1:3), '../tex/tables/bootstrap1_1.tex')
drosPrintBootstrapResults(res1(:, :, 4:6), 1:4, {'MAP', 'ML', 'Reg', 'Inf'}, T1(4:6), '../tex/tables/bootstrap1_2.tex')

drosPrintBootstrapResults(res1c(:, :, 1:3), 1:4, {'MAP', 'ML', 'Reg', 'Inf'}, T1(1:3), '../tex/tables/bootstrap1c_1.tex')
drosPrintBootstrapResults(res1c(:, :, 4:6), 1:4, {'MAP', 'ML', 'Reg', 'Inf'}, T1(4:6), '../tex/tables/bootstrap1c_2.tex')

drosPrintBootstrapResults(res2(:, :, 1:3), 1:4, {'P32', 'P2', 'ML', 'Inf'}, T2(1:3), '../tex/tables/bootstrap2_1.tex')
drosPrintBootstrapResults(res2(:, :, 4:6), 1:4, {'P32', 'P2', 'ML', 'Inf'}, T2(4:6), '../tex/tables/bootstrap2_2.tex')
drosPrintBootstrapResults(res2(:, :, 5:7), 1:4, {'P32', 'P2', 'ML', 'Inf'}, T2(5:7), '../tex/tables/bootstrap2_3.tex')

drosPrintBootstrapResults(res3(:, :, 1:3), 1:4, {'P32', 'P2', 'ML', 'Inf'}, T2(1:3), '../tex/tables/bootstrap3_1.tex')
drosPrintBootstrapResults(res3(:, :, 4:6), 1:4, {'P32', 'P2', 'ML', 'Inf'}, T2(4:6), '../tex/tables/bootstrap3_2.tex')
drosPrintBootstrapResults(res3(:, :, 5:7), 1:4, {'P32', 'P2', 'ML', 'Inf'}, T2(5:7), '../tex/tables/bootstrap3_3.tex')
