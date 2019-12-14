module choice_data
  !
  ! v1.0 17.04.01  Date written.
  !    
  ! Contains default data for components and tabulated data for noise empirical models.
  !
  use sim_precision
  use units_and_constants
  !
  implicit none
  !
  real(kind=rp), parameter :: dTref = 0.5550_rp
  real(kind=rp), parameter ::  mref = 0.4536_rp
  !
  real(kind=rp), parameter :: p0 = 2.0000E-05_rp
  real(kind=rp), parameter :: Tisa = 288.15_rp ! ISA Condition Static Temperature under ISA, sea level condition
  real(kind=rp), parameter :: Pisa = 101325.0_rp ! ISA Condition Static Pressure under ISA, sea level condition
  real(kind=rp), parameter :: Risa = 287.05_rp ! ISA Condition gas constant
  real(kind=rp), parameter :: rhoisa = Pisa/(Risa*Tisa) ! ISA Condition Static Density under ISA, sea level condition  
  !
  integer, parameter :: n_dir_data = 37 ! possible generalize all
  integer, parameter :: n_rows = 36
  integer, parameter :: n_cols = 7
  real(kind=rp), dimension(n_rows,n_cols+1), parameter :: table_iv = (/ & 
         & -41.64_rp,64.84_rp,5.85_rp,45.66_rp,-45.95_rp,28.98_rp,-67.28_rp,65.60_rp,-3.78_rp, & 
         & -3.81_rp,-12.14_rp,-6.75_rp,14.77_rp,37.95_rp,0.75_rp,-11.83_rp,-6.53_rp,12.77_rp, & 
         & -17.63_rp,5.05_rp,-0.12_rp,11.86_rp,-41.33_rp,40.51_rp,-37.63_rp,29.26_rp,33.74_rp, & 
         & -12.19_rp,68.07_rp,-69.25_rp,23.03_rp,10.90_rp,-1.15_rp,34.02_rp,8.42_rp,43.73_rp, & 
         & -12.43_rp,-17.35_rp,3.35_rp,5.25_rp,-75.93_rp,142.80_rp,-147.56_rp,180.53_rp,10.22_rp, & 
         & -2.33_rp,9.26_rp,15.17_rp,-11.61_rp,0.47_rp,-53.84_rp,3.31_rp,43.94_rp,-43.15_rp, & 
         & 3.14_rp,44.33_rp,-4.68_rp,37.70_rp,-1.70_rp,18.73_rp,-23.67_rp,-0.44_rp,-43.67_rp, & 
         & -13.22_rp,19.20_rp,-7.82_rp,-17.62_rp,50.11_rp,23.04_rp,-2.11_rp,3.65_rp,-8.43_rp, & 
         & -10.57_rp,-16.34_rp,0.54_rp,-6.42_rp,-56.72_rp,107.84_rp,-178.07_rp,168.22_rp,13.12_rp, & 
         & -3.30_rp,8.54_rp,14.35_rp,-1.72_rp,-1.10_rp,-52.80_rp,3.71_rp,46.06_rp,-34.62_rp, & 
         & 3.37_rp,25.17_rp,-2.74_rp,13.31_rp,-19.93_rp,18.16_rp,-28.99_rp,-2.44_rp,-32.07_rp, & 
         & -27.97_rp,20.33_rp,0.94_rp,-12.62_rp,32.54_rp,11.21_rp,-4.00_rp,19.64_rp,-6.21_rp, & 
         & -7.84_rp,-16.71_rp,1.26_rp,11.35_rp,-47.74_rp,76.07_rp,-192.58_rp,119.55_rp,13.54_rp, & 
         & -3.39_rp,3.45_rp,5.40_rp,13.00_rp,0.91_rp,-39.79_rp,1.46_rp,42.20_rp,-24.48_rp, & 
         & 4.10_rp,6.61_rp,-4.62_rp,5.29_rp,-59.88_rp,12.75_rp,-27.16_rp,0.44_rp,-13.20_rp, & 
         & -54.74_rp,23.72_rp,-1.97_rp,-5.29_rp,26.46_rp,4.20_rp,-1.86_rp,27.43_rp,7.00_rp, & 
         & -4.54_rp,-12.48_rp,0.20_rp,-28.69_rp,-14.32_rp,28.13_rp,-137.39_rp,43.89_rp,9.91_rp, & 
         & -0.83_rp,-2.45_rp,15.72_rp,8.69_rp,3.07_rp,-14.75_rp,1.27_rp,33.04_rp,-11.60_rp, & 
         & 9.30_rp,-1.50_rp,-2.88_rp,-17.63_rp,-16.95_rp,9.76_rp,-35.00_rp,-1.72_rp,-29.24_rp, & 
         & -37.21_rp,23.68_rp,-4.61_rp,3.08_rp,15.57_rp,-4.62_rp,-6.54_rp,10.65_rp,-2.40_rp, & 
         & 0.44_rp,-4.83_rp,2.24_rp,-21.59_rp,2.43_rp,13.19_rp,-83.21_rp,55.58_rp,5.39_rp, & 
         & -1.46_rp,-1.91_rp,10.19_rp,9.09_rp,5.28_rp,-20.67_rp,2.11_rp,16.57_rp,0.30_rp, & 
         & 5.16_rp,1.27_rp,-11.52_rp,-30.21_rp,-22.64_rp,2.25_rp,-11.17_rp,-8.14_rp,-21.80_rp, & 
         & -48.71_rp,11.00_rp,6.05_rp,1.40_rp,26.99_rp,9.05_rp,-0.49_rp,3.77_rp,21.33_rp, & 
         & 6.29_rp,6.87_rp,-0.20_rp,-19.65_rp,8.96_rp,-24.37_rp,21.51_rp,-30.98_rp,-4.83_rp, & 
         & -0.42_rp,2.45_rp,-1.38_rp,-14.63_rp,-4.27_rp,17.42_rp,-1.23_rp,-18.83_rp,12.34_rp, & 
         & -4.12_rp,-14.60_rp,9.39_rp,5.37_rp,30.05_rp,7.51_rp,-2.11_rp,7.46_rp,6.99_rp, & 
         & 52.23_rp,-20.08_rp,8.12_rp,-0.34_rp,-8.18_rp,-14.19_rp,-3.57_rp,10.43_rp,-14.97_rp, & 
         & 6.38_rp,4.36_rp,-12.57_rp,-51.73_rp,-14.87_rp,-69.90_rp,102.77_rp,179.12_rp,-12.40_rp, & 
         & 9.62_rp,60.26_rp,9.77_rp,-1.29_rp,-19.71_rp,56.67_rp,-19.54_rp,-24.69_rp,16.01_rp, & 
         & 31.81_rp,-44.19_rp,-5.09_rp,148.29_rp,111.51_rp,105.11_rp,-14.68_rp,69.81_rp,52.26_rp, & 
         & 111.33_rp,-13.36_rp,-8.01_rp,49.66_rp,-20.66_rp,-60.82_rp,-25.51_rp,-44.59_rp,-11.41_rp/)
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_v = (/ & 
         & -27.59_rp,-8.46_rp,19.20_rp,0.00_rp,0.00_rp,-40.02_rp,0.00_rp,0.00_rp,13.65_rp, & 
         & -10.32_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,19.06_rp, & 
         & 4.30_rp,-9.09_rp,-36.86_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp, & 
         & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp, & 
         & -14.16_rp,8.44_rp,-0.12_rp,25.92_rp,-47.71_rp,15.36_rp,-138.82_rp,103.69_rp,2.61_rp, & 
         & 2.16_rp,10.82_rp,8.90_rp,-19.14_rp,-15.81_rp,-6.02_rp,2.97_rp,17.79_rp,-11.70_rp, & 
         & 42.09_rp,-27.74_rp,-5.00_rp,0.00_rp,2.81_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp, & 
         & 0.00_rp,-72.96_rp,98.39_rp,-18.80_rp,-74.92_rp,-6.34_rp,-75.38_rp,13.14_rp,-1.47_rp, & 
         & -10.17_rp,1.49_rp,-0.77_rp,-3.57_rp,-7.41_rp,-3.01_rp,56.76_rp,46.69_rp,-12.30_rp, & 
         & 1.81_rp,-24.04_rp,12.25_rp,-4.07_rp,-17.73_rp,19.18_rp,0.77_rp,-35.24_rp,15.35_rp, & 
         & 1.49_rp,-21.10_rp,6.32_rp,-14.08_rp,57.68_rp,-48.21_rp,9.70_rp,-8.41_rp,-6.49_rp, & 
         & 27.94_rp,-64.83_rp,52.29_rp,12.53_rp,-8.12_rp,-41.86_rp,-5.52_rp,25.61_rp,-13.69_rp, & 
         & -13.97_rp,-4.45_rp,-1.52_rp,-4.99_rp,27.91_rp,-16.45_rp,113.51_rp,-28.24_rp,-1.04_rp, & 
         & 1.38_rp,-5.18_rp,6.03_rp,5.56_rp,-10.58_rp,12.39_rp,-1.34_rp,-10.08_rp,15.94_rp, & 
         & -22.60_rp,-9.46_rp,0.04_rp,-30.83_rp,11.72_rp,-14.29_rp,39.73_rp,-8.39_rp,5.86_rp, & 
         & -6.40_rp,7.96_rp,5.73_rp,4.15_rp,-8.41_rp,-5.54_rp,32.27_rp,16.43_rp,12.21_rp, & 
         & -17.56_rp,-6.33_rp,-8.44_rp,2.30_rp,21.08_rp,-0.99_rp,148.84_rp,60.35_rp,11.00_rp, & 
         & 3.47_rp,38.23_rp,-12.64_rp,8.09_rp,-19.32_rp,5.69_rp,2.86_rp,29.79_rp,-8.78_rp, & 
         & -28.44_rp,24.30_rp,3.37_rp,11.15_rp,14.29_rp,40.54_rp,63.43_rp,-1.51_rp,-1.82_rp, & 
         & -35.42_rp,81.68_rp,-32.15_rp,12.02_rp,-42.27_rp,15.64_rp,62.84_rp,-34.37_rp,-18.90_rp, & 
         & -23.89_rp,-20.25_rp,-1.12_rp,-24.48_rp,-6.53_rp,81.36_rp,-87.26_rp,308.24_rp,15.90_rp, & 
         & 1.82_rp,5.54_rp,-22.06_rp,4.88_rp,-6.16_rp,18.48_rp,6.38_rp,83.56_rp,-39.83_rp, & 
         & -9.27_rp,36.32_rp,3.23_rp,-67.15_rp,-48.23_rp,72.22_rp,-13.79_rp,5.93_rp,-40.51_rp, & 
         & -78.46_rp,116.27_rp,-35.47_rp,-3.24_rp,-33.64_rp,42.59_rp,44.41_rp,-27.71_rp,-17.02_rp, & 
         & -29.75_rp,-44.32_rp,-20.97_rp,0.00_rp,-51.27_rp,0.00_rp,0.00_rp,0.00_rp,10.24_rp, & 
         & 2.51_rp,-25.51_rp,0.00_rp,60.14_rp,0.00_rp,0.00_rp,-1.02_rp,189.29_rp,-51.46_rp, & 
         & 25.99_rp,17.11_rp,35.87_rp,0.00_rp,0.00_rp,83.44_rp,0.00_rp,57.28_rp,0.00_rp, & 
         & -57.85_rp,187.73_rp,0.00_rp,-42.41_rp,0.00_rp,0.00_rp,-51.21_rp,211.60_rp,0.00_rp /)
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_1 = (/ & 
         & 3.70_rp,7.34_rp,-10.40_rp,0.00_rp,0.00_rp,-18.95_rp,0.00_rp,0.00_rp,-14.34_rp,   & 
         & 7.58_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,-56.30_rp,   & 
         & 23.01_rp,116.22_rp,-153.51_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & -2.47_rp,-23.49_rp,12.67_rp,-4.52_rp,79.23_rp,-80.04_rp,129.39_rp,-110.57_rp,5.46_rp,   & 
         & -8.73_rp,-1.96_rp,-41.58_rp,-3.23_rp,0.99_rp,4.78_rp,11.20_rp,-30.45_rp,72.74_rp,   & 
         & -20.76_rp,-68.87_rp,-14.73_rp,0.00_rp,-52.60_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,-40.89_rp,28.94_rp,-14.94_rp,97.14_rp,-28.07_rp,-1.42_rp,52.79_rp,16.83_rp,   & 
         & -1.86_rp,1.68_rp,3.69_rp,45.78_rp,1.83_rp,-26.44_rp,-56.85_rp,-133.33_rp,11.36_rp,   & 
         & -6.40_rp,-6.81_rp,29.32_rp,61.90_rp,46.61_rp,-90.89_rp,-14.99_rp,13.60_rp,29.52_rp,   & 
         & 28.84_rp,-49.73_rp,-52.93_rp,-5.21_rp,-32.86_rp,-24.66_rp,30.99_rp,-0.45_rp,8.72_rp,   & 
         & -112.81_rp,3.99_rp,9.23_rp,-6.12_rp,89.70_rp,27.76_rp,-37.32_rp,43.89_rp,114.16_rp,   & 
         & 2.47_rp,3.78_rp,-2.71_rp,-58.04_rp,13.59_rp,-13.19_rp,-13.42_rp,50.04_rp,-7.62_rp,   & 
         & 7.93_rp,4.70_rp,13.68_rp,-17.86_rp,-18.29_rp,17.67_rp,7.15_rp,-24.54_rp,2.09_rp,   & 
         & 3.56_rp,-0.92_rp,5.43_rp,-7.88_rp,34.90_rp,1.03_rp,4.62_rp,9.05_rp,3.72_rp,   & 
         & 54.06_rp,-18.37_rp,-6.37_rp,16.18_rp,-30.99_rp,-38.71_rp,-3.53_rp,-10.51_rp,-27.20_rp,   & 
         & 3.43_rp,-7.11_rp,-1.73_rp,-134.84_rp,30.92_rp,-11.66_rp,-148.64_rp,517.96_rp,-27.22_rp,   & 
         & 13.60_rp,-10.47_rp,-21.77_rp,-75.31_rp,-68.09_rp,100.51_rp,22.38_rp,-57.11_rp,-20.12_rp,   & 
         & -37.31_rp,34.58_rp,65.15_rp,-27.75_rp,78.78_rp,11.79_rp,-51.20_rp,5.80_rp,1.63_rp,   & 
         & 155.31_rp,43.95_rp,11.03_rp,16.62_rp,-50.66_rp,-55.94_rp,54.04_rp,-3.60_rp,-133.30_rp,   & 
         & 3.09_rp,-3.09_rp,7.00_rp,-144.66_rp,19.66_rp,48.07_rp,-310.26_rp,444.11_rp,-35.73_rp,   & 
         & 15.86_rp,-32.32_rp,-34.43_rp,-82.60_rp,-46.64_rp,145.35_rp,19.47_rp,-93.90_rp,-23.78_rp,   & 
         & -39.03_rp,75.72_rp,87.70_rp,-91.66_rp,6.24_rp,-24.29_rp,-101.39_rp,-38.76_rp,42.76_rp,   & 
         & 144.63_rp,-41.88_rp,-57.54_rp,25.73_rp,-67.03_rp,-78.03_rp,26.98_rp,-27.94_rp,-135.23_rp,   & 
         & 0.79_rp,-5.46_rp,10.52_rp,0.00_rp,145.06_rp,0.00_rp,0.00_rp,0.00_rp,-16.11_rp,   & 
         & 9.31_rp,21.79_rp,0.00_rp,55.93_rp,0.00_rp,0.00_rp,35.58_rp,46.63_rp,-42.10_rp,   & 
         & 17.53_rp,6.62_rp,24.43_rp,0.00_rp,0.00_rp,11.27_rp,0.00_rp,-27.01_rp,0.00_rp,   & 
         & -58.70_rp,118.17_rp,0.00_rp,40.40_rp,0.00_rp,0.00_rp,-10.95_rp,253.16_rp,0.00_rp/) 
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_2 = (/ & 
         & 1.22_rp,3.17_rp,-0.48_rp,0.00_rp,0.00_rp,-26.42_rp,0.00_rp,0.00_rp,-9.00_rp,   & 
         & 3.40_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,-40.64_rp,   & 
         & 21.42_rp,110.97_rp,-134.11_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & -3.65_rp,-12.44_rp,11.92_rp,16.82_rp,29.00_rp,-24.89_rp,-65.02_rp,-61.09_rp,-1.28_rp,   & 
         & -3.79_rp,-9.73_rp,-35.24_rp,-23.90_rp,-10.10_rp,26.77_rp,16.09_rp,-37.30_rp,52.12_rp,   & 
         & -32.66_rp,-36.49_rp,-9.04_rp,0.00_rp,40.74_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,-10.65_rp,-6.73_rp,7.83_rp,79.16_rp,-33.59_rp,34.85_rp,29.77_rp,38.85_rp,   & 
         & -2.37_rp,1.37_rp,5.55_rp,31.26_rp,-6.71_rp,-1.74_rp,-98.78_rp,-90.94_rp,9.15_rp,   & 
         & -5.02_rp,-9.05_rp,30.51_rp,42.28_rp,43.55_rp,-77.65_rp,-10.81_rp,9.24_rp,24.29_rp,   & 
         & 25.96_rp,-37.40_rp,-45.04_rp,-15.31_rp,-35.61_rp,-27.28_rp,32.81_rp,-2.58_rp,3.39_rp,   & 
         & -82.89_rp,0.48_rp,-15.65_rp,-9.92_rp,82.58_rp,27.01_rp,-37.25_rp,28.51_rp,93.89_rp,   & 
         & 2.62_rp,2.62_rp,-1.35_rp,-26.27_rp,-2.08_rp,-8.55_rp,-21.59_rp,38.34_rp,-6.66_rp,   & 
         & 6.10_rp,5.80_rp,5.86_rp,-9.11_rp,-14.30_rp,11.22_rp,2.55_rp,-23.27_rp,2.35_rp,   & 
         & 4.99_rp,-5.28_rp,2.02_rp,15.41_rp,16.35_rp,1.25_rp,7.86_rp,15.43_rp,14.65_rp,   & 
         & 44.36_rp,-20.26_rp,-7.14_rp,15.08_rp,-16.82_rp,-35.62_rp,-4.91_rp,-0.75_rp,-14.88_rp,   & 
         & 4.47_rp,-2.10_rp,-4.34_rp,-60.59_rp,21.60_rp,-50.63_rp,73.99_rp,266.55_rp,-23.97_rp,   & 
         & 7.71_rp,11.00_rp,-25.82_rp,-58.50_rp,-60.97_rp,54.51_rp,18.98_rp,-57.94_rp,-6.73_rp,   & 
         & -27.44_rp,8.84_rp,36.34_rp,45.05_rp,66.24_rp,22.95_rp,-11.47_rp,24.57_rp,12.79_rp,   & 
         & 134.04_rp,-60.13_rp,20.79_rp,9.30_rp,19.56_rp,-40.64_rp,41.64_rp,16.14_rp,-82.56_rp,   & 
         & 4.41_rp,0.79_rp,4.23_rp,-78.79_rp,29.71_rp,-18.07_rp,-88.84_rp,177.25_rp,-30.07_rp,   & 
         & 8.90_rp,-14.57_rp,-35.86_rp,-47.12_rp,-41.78_rp,92.39_rp,9.40_rp,-84.83_rp,-7.98_rp,   & 
         & -29.21_rp,41.15_rp,53.58_rp,-46.19_rp,1.65_rp,-11.08_rp,-62.92_rp,-17.51_rp,52.01_rp,   & 
         & 103.78_rp,-56.46_rp,-30.57_rp,14.24_rp,0.54_rp,-55.15_rp,15.97_rp,-2.06_rp,-73.93_rp,   & 
         & 1.94_rp,-3.77_rp,1.59_rp,0.00_rp,48.83_rp,0.00_rp,0.00_rp,0.00_rp,-14.78_rp,   & 
         & 4.99_rp,7.99_rp,0.00_rp,25.90_rp,0.00_rp,0.00_rp,28.87_rp,18.83_rp,-35.89_rp,   & 
         & -0.04_rp,1.35_rp,8.25_rp,0.00_rp,0.00_rp,26.79_rp,0.00_rp,-3.62_rp,0.00_rp,   & 
         & -25.88_rp,37.24_rp,0.00_rp,7.52_rp,0.00_rp,0.00_rp,3.23_rp,140.03_rp,0.00_rp/)
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_3 = (/ & 
         & -1.04_rp,1.68_rp,8.40_rp,0.00_rp,0.00_rp,-54.75_rp,0.00_rp,0.00_rp,-8.24_rp,   & 
         & 1.85_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,-34.67_rp,   & 
         & 15.08_rp,90.30_rp,-94.68_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & -4.01_rp,-10.40_rp,8.03_rp,-15.22_rp,32.16_rp,-22.01_rp,67.03_rp,-63.85_rp,-3.07_rp,   & 
         & -4.67_rp,-8.12_rp,-21.02_rp,-17.50_rp,-6.67_rp,3.70_rp,11.84_rp,-34.77_rp,43.27_rp,   & 
         & -27.72_rp,-25.73_rp,-6.51_rp,0.00_rp,84.30_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,-28.94_rp,8.84_rp,-5.87_rp,94.57_rp,-18.78_rp,25.00_rp,22.52_rp,15.45_rp,   & 
         & -2.59_rp,0.82_rp,5.49_rp,15.47_rp,-10.11_rp,11.27_rp,-98.32_rp,-73.38_rp,7.45_rp,   & 
         & -3.65_rp,-10.09_rp,26.33_rp,23.21_rp,33.87_rp,-56.95_rp,-7.34_rp,9.51_rp,12.71_rp,   & 
         & 19.30_rp,-22.37_rp,-28.63_rp,-19.09_rp,-17.47_rp,-24.30_rp,21.20_rp,-5.30_rp,-2.59_rp,   & 
         & -46.18_rp,-0.50_rp,-13.79_rp,-12.63_rp,56.22_rp,24.01_rp,-31.07_rp,15.73_rp,59.52_rp,   & 
         & 2.49_rp,3.27_rp,-0.05_rp,-10.56_rp,-10.91_rp,-3.47_rp,-37.55_rp,15.96_rp,-7.22_rp,   & 
         & 5.41_rp,1.40_rp,3.39_rp,-13.25_rp,-14.76_rp,16.79_rp,2.45_rp,-22.56_rp,-0.88_rp,   & 
         & 4.08_rp,-3.00_rp,6.53_rp,19.66_rp,-0.33_rp,-2.10_rp,-0.69_rp,15.67_rp,14.39_rp,   & 
         & 47.14_rp,-25.94_rp,-0.10_rp,12.37_rp,-20.01_rp,-33.64_rp,-5.05_rp,-0.36_rp,-20.12_rp,   & 
         & 4.71_rp,2.65_rp,-2.26_rp,-16.55_rp,4.78_rp,-40.41_rp,78.35_rp,127.54_rp,-21.48_rp,   & 
         & 5.29_rp,8.36_rp,-20.13_rp,-46.79_rp,-48.93_rp,36.10_rp,16.31_rp,-54.45_rp,-1.44_rp,   & 
         & -20.55_rp,-0.65_rp,21.71_rp,52.86_rp,18.33_rp,19.39_rp,-7.37_rp,29.17_rp,10.40_rp,   & 
         & 97.73_rp,-62.34_rp,23.46_rp,7.23_rp,38.75_rp,-31.31_rp,33.85_rp,18.51_rp,-52.72_rp,   & 
         & 5.31_rp,4.65_rp,0.86_rp,-10.41_rp,24.80_rp,-60.99_rp,104.97_rp,-46.86_rp,-21.52_rp,   & 
         & 3.11_rp,3.94_rp,-32.52_rp,-11.60_rp,-29.83_rp,50.35_rp,0.34_rp,-71.08_rp,8.36_rp,   & 
         & -20.51_rp,10.03_rp,23.89_rp,14.42_rp,18.42_rp,-7.09_rp,-13.42_rp,-2.62_rp,55.69_rp,   & 
         & 59.63_rp,-56.14_rp,-14.60_rp,7.76_rp,33.84_rp,-36.27_rp,9.21_rp,8.29_rp,-22.75_rp,   & 
         & 5.70_rp,13.61_rp,-1.43_rp,0.00_rp,-18.04_rp,0.00_rp,0.00_rp,0.00_rp,-3.63_rp,   & 
         & -5.95_rp,15.90_rp,0.00_rp,-13.52_rp,0.00_rp,0.00_rp,22.72_rp,-30.83_rp,-10.62_rp,   & 
         & -27.21_rp,-9.68_rp,-5.39_rp,0.00_rp,0.00_rp,8.13_rp,0.00_rp,-23.46_rp,0.00_rp,   & 
         & -7.09_rp,-30.64_rp,0.00_rp,-3.71_rp,0.00_rp,0.00_rp,24.65_rp,-23.59_rp,0.00_rp  /) 
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_4 = (/ & 
         & -1.64_rp,0.01_rp,-10.22_rp,0.00_rp,0.00_rp,55.56_rp,0.00_rp,0.00_rp,-9.26_rp,   & 
         & 2.52_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,-24.61_rp,   & 
         & 2.26_rp,30.25_rp,-13.68_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & -5.07_rp,-5.84_rp,2.57_rp,-3.41_rp,24.45_rp,-4.72_rp,121.57_rp,-115.10_rp,-7.37_rp,   & 
         & -2.83_rp,-6.49_rp,-7.52_rp,0.15_rp,-6.25_rp,-9.20_rp,6.09_rp,-34.43_rp,14.69_rp,   & 
         & -29.04_rp,7.18_rp,7.44_rp,0.00_rp,52.85_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,-14.76_rp,-27.12_rp,-7.97_rp,97.36_rp,-5.12_rp,32.30_rp,9.93_rp,-17.83_rp,   & 
         & -2.41_rp,-1.71_rp,4.45_rp,-0.09_rp,-5.45_rp,9.32_rp,-71.18_rp,-44.51_rp,4.64_rp,   & 
         & -2.47_rp,-3.32_rp,9.00_rp,2.08_rp,11.58_rp,-28.32_rp,-3.78_rp,6.60_rp,-1.76_rp,   & 
         & 8.14_rp,-8.14_rp,-9.72_rp,-20.17_rp,-9.48_rp,-16.28_rp,14.87_rp,-6.62_rp,5.94_rp,   & 
         & -0.69_rp,-9.06_rp,-3.05_rp,-16.56_rp,29.94_rp,19.28_rp,-20.86_rp,8.12_rp,26.20_rp,   & 
         & 2.46_rp,4.05_rp,1.99_rp,-5.26_rp,-17.12_rp,1.67_rp,-90.44_rp,-4.80_rp,-5.91_rp,   & 
         & 2.04_rp,-2.26_rp,3.32_rp,-13.37_rp,-10.70_rp,10.92_rp,1.41_rp,-17.84_rp,-0.14_rp,   & 
         & 9.82_rp,-9.00_rp,4.43_rp,9.20_rp,-16.29_rp,-5.91_rp,-15.15_rp,11.31_rp,9.25_rp,   & 
         & 38.42_rp,-40.22_rp,21.81_rp,7.49_rp,-12.12_rp,-27.67_rp,-16.58_rp,8.77_rp,-9.70_rp,   & 
         & 4.07_rp,8.19_rp,2.55_rp,4.20_rp,-3.89_rp,-11.66_rp,-28.79_rp,-9.91_rp,-11.96_rp,   & 
         & 1.12_rp,-11.98_rp,-1.24_rp,-14.47_rp,-12.29_rp,16.60_rp,-1.65_rp,-29.55_rp,10.36_rp,   & 
         & 5.10_rp,-17.37_rp,0.43_rp,7.91_rp,-19.67_rp,-1.19_rp,-35.69_rp,16.32_rp,1.63_rp,   & 
         & 28.92_rp,-56.62_rp,46.23_rp,7.59_rp,6.42_rp,-19.83_rp,-5.19_rp,12.36_rp,1.61_rp,   & 
         & 5.42_rp,7.07_rp,0.41_rp,-16.48_rp,6.16_rp,-35.36_rp,20.68_rp,-70.58_rp,-8.03_rp,   & 
         & 0.34_rp,6.59_rp,-6.76_rp,11.29_rp,-2.58_rp,30.37_rp,-4.28_rp,-49.04_rp,23.52_rp,   & 
         & -13.15_rp,-10.68_rp,10.26_rp,17.34_rp,57.36_rp,-30.38_rp,18.35_rp,-13.90_rp,30.34_rp,   & 
         & 19.82_rp,-34.81_rp,-11.04_rp,9.20_rp,-9.61_rp,-30.84_rp,2.62_rp,-0.88_rp,-6.07_rp,   & 
         & 5.12_rp,19.60_rp,8.46_rp,0.00_rp,29.34_rp,0.00_rp,0.00_rp,0.00_rp,8.45_rp,   & 
         & -6.16_rp,23.63_rp,0.00_rp,-8.09_rp,0.00_rp,0.00_rp,-7.00_rp,-55.14_rp,28.75_rp,   & 
         & -35.69_rp,3.12_rp,-0.50_rp,0.00_rp,0.00_rp,-48.93_rp,0.00_rp,-56.75_rp,0.00_rp,   & 
         & -6.66_rp,-8.59_rp,0.00_rp,22.75_rp,0.00_rp,0.00_rp,33.55_rp,-64.11_rp,0.00_rp/)
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_5 = (/ & 
         & -3.75_rp,-10.99_rp,6.77_rp,0.00_rp,0.00_rp,-27.09_rp,0.00_rp,0.00_rp,-1.58_rp,   & 
         & -0.66_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,3.96_rp,   & 
         & -3.56_rp,-28.10_rp,28.47_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & -5.07_rp,-5.83_rp,6.36_rp,-1.46_rp,8.00_rp,20.42_rp,-17.17_rp,-82.36_rp,-4.94_rp,   & 
         & -3.83_rp,-0.15_rp,9.98_rp,18.91_rp,-0.31_rp,-9.64_rp,-7.50_rp,-23.61_rp,1.57_rp,   & 
         & -14.05_rp,6.69_rp,3.23_rp,0.00_rp,-13.92_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
         & 0.00_rp,2.05_rp,-62.76_rp,-10.84_rp,97.85_rp,-2.13_rp,-1.19_rp,-2.44_rp,11.61_rp,   & 
         & -1.38_rp,-0.71_rp,4.60_rp,5.42_rp,-5.25_rp,4.04_rp,-79.12_rp,-4.86_rp,0.77_rp,   & 
         & -1.68_rp,-16.01_rp,9.60_rp,-1.24_rp,2.47_rp,-9.53_rp,4.34_rp,3.52_rp,-2.95_rp,   & 
         & 6.19_rp,-6.91_rp,-0.63_rp,-39.90_rp,-18.90_rp,-24.09_rp,-1.02_rp,-15.83_rp,-8.10_rp,   & 
         & -25.06_rp,-5.22_rp,-5.59_rp,-8.15_rp,43.33_rp,8.79_rp,-14.11_rp,14.15_rp,10.29_rp,   & 
         & 2.84_rp,5.74_rp,-0.15_rp,-1.37_rp,-23.25_rp,-3.05_rp,-92.67_rp,-1.39_rp,-3.71_rp,   & 
         & -0.59_rp,1.83_rp,0.45_rp,-9.72_rp,-4.45_rp,-1.42_rp,2.51_rp,-12.23_rp,3.88_rp,   & 
         & 14.59_rp,-7.38_rp,0.65_rp,14.80_rp,2.66_rp,-4.64_rp,-16.51_rp,4.07_rp,4.35_rp,   & 
         & 28.55_rp,-37.38_rp,23.58_rp,5.74_rp,10.87_rp,-20.43_rp,-24.01_rp,3.45_rp,-6.60_rp,   & 
         & 2.78_rp,8.91_rp,-1.05_rp,9.15_rp,-0.79_rp,-23.47_rp,13.61_rp,-116.57_rp,-4.01_rp,   & 
         & -2.15_rp,-7.76_rp,-2.78_rp,3.36_rp,9.74_rp,4.99_rp,-7.77_rp,-12.54_rp,19.95_rp,   & 
         & 16.09_rp,-26.65_rp,-3.58_rp,35.21_rp,16.38_rp,5.10_rp,-54.63_rp,16.82_rp,13.98_rp,   & 
         & 31.90_rp,-47.57_rp,54.21_rp,6.65_rp,-15.19_rp,-18.55_rp,-25.52_rp,9.84_rp,10.36_rp,   & 
         & 2.90_rp,8.99_rp,-5.86_rp,9.47_rp,23.37_rp,-60.23_rp,219.91_rp,-264.53_rp,3.44_rp,   & 
         & -3.59_rp,23.14_rp,-4.94_rp,14.36_rp,11.03_rp,14.15_rp,-12.63_rp,-20.67_rp,25.73_rp,   & 
         & -5.74_rp,-26.72_rp,7.60_rp,74.22_rp,80.70_rp,-2.84_rp,23.63_rp,6.56_rp,33.27_rp,   & 
         & 45.37_rp,-7.72_rp,-7.63_rp,10.65_rp,-46.29_rp,-24.10_rp,-0.60_rp,-3.12_rp,-6.70_rp,   & 
         & 2.00_rp,22.10_rp,8.12_rp,0.00_rp,74.06_rp,0.00_rp,0.00_rp,0.00_rp,9.85_rp,   & 
         & -7.35_rp,30.71_rp,0.00_rp,-9.17_rp,0.00_rp,0.00_rp,-12.49_rp,-34.22_rp,32.38_rp,   & 
         & -12.25_rp,-20.14_rp,-2.19_rp,0.00_rp,0.00_rp,-44.56_rp,0.00_rp,-50.69_rp,0.00_rp,   & 
         & 20.46_rp,32.20_rp,0.00_rp,35.61_rp,0.00_rp,0.00_rp,6.02_rp,-40.05_rp,0.00_rp /)
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_6 = (/ & 
   & 0.57_rp,-3.13_rp,-7.24_rp,0.00_rp,0.00_rp,68.97_rp,0.00_rp,0.00_rp,1.92_rp,   & 
   & 2.55_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,7.66_rp,   & 
   & 7.72_rp,1.75_rp,-33.12_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
   & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
   & 1.48_rp,-3.25_rp,3.72_rp,-20.25_rp,11.92_rp,-4.15_rp,12.63_rp,-67.21_rp,1.97_rp,   & 
   & -2.33_rp,2.84_rp,-4.78_rp,15.12_rp,0.80_rp,-11.26_rp,0.62_rp,-3.72_rp,-0.22_rp,   & 
   & -10.19_rp,14.22_rp,0.31_rp,0.00_rp,-51.56_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
   & 0.00_rp,21.93_rp,-53.50_rp,-1.44_rp,79.25_rp,-1.08_rp,7.52_rp,3.99_rp,12.71_rp,   & 
   & 1.16_rp,-1.81_rp,-0.78_rp,-1.80_rp,11.21_rp,-22.85_rp,31.88_rp,-44.76_rp,1.92_rp,   & 
   & -1.43_rp,-2.80_rp,3.29_rp,2.02_rp,-4.05_rp,-9.27_rp,2.58_rp,11.77_rp,-4.35_rp,   & 
   & 0.37_rp,2.96_rp,-1.74_rp,23.56_rp,4.43_rp,2.85_rp,-25.37_rp,7.39_rp,-24.99_rp,   & 
   & -6.14_rp,3.32_rp,17.20_rp,0.70_rp,0.40_rp,6.44_rp,9.15_rp,-6.95_rp,-10.21_rp,   & 
   & -4.05_rp,3.84_rp,-9.83_rp,43.37_rp,41.01_rp,-30.81_rp,242.06_rp,-99.58_rp,0.41_rp,   & 
   & -0.21_rp,16.32_rp,-5.57_rp,31.06_rp,-2.26_rp,7.92_rp,-13.42_rp,-3.41_rp,11.95_rp,   & 
   & -6.48_rp,9.16_rp,-3.54_rp,14.44_rp,27.26_rp,7.28_rp,38.64_rp,-0.14_rp,37.73_rp,   & 
   & -24.67_rp,35.00_rp,-31.01_rp,17.49_rp,11.73_rp,-5.40_rp,16.30_rp,-32.98_rp,9.01_rp,   & 
   & -8.48_rp,11.45_rp,-7.76_rp,69.04_rp,40.45_rp,-25.27_rp,152.09_rp,-224.10_rp,23.11_rp,   & 
   & -6.60_rp,-10.12_rp,10.57_rp,34.53_rp,53.96_rp,-21.05_rp,-25.25_rp,56.88_rp,9.86_rp,   & 
   & 13.93_rp,-6.88_rp,-12.28_rp,-3.67_rp,67.90_rp,-22.91_rp,15.03_rp,-25.53_rp,2.31_rp,   & 
   & -30.80_rp,93.37_rp,-67.11_rp,-3.57_rp,-9.13_rp,23.62_rp,-24.65_rp,-36.57_rp,29.39_rp,   & 
   & -9.64_rp,8.21_rp,-3.31_rp,48.36_rp,46.86_rp,18.50_rp,115.90_rp,-68.30_rp,16.02_rp,   & 
   & -0.98_rp,-18.37_rp,7.82_rp,10.95_rp,20.53_rp,7.66_rp,-22.41_rp,22.00_rp,19.70_rp,   & 
   & 5.61_rp,0.10_rp,-5.68_rp,-79.52_rp,41.75_rp,-68.32_rp,103.03_rp,-26.80_rp,26.87_rp,   & 
   & 21.39_rp,55.56_rp,-66.72_rp,3.56_rp,-13.06_rp,-13.87_rp,-26.73_rp,-45.42_rp,17.86_rp,   & 
   & -6.38_rp,-0.40_rp,10.60_rp,0.00_rp,45.50_rp,0.00_rp,0.00_rp,0.00_rp,-2.13_rp,   & 
   & -1.19_rp,-19.07_rp,0.00_rp,10.60_rp,0.00_rp,0.00_rp,-5.96_rp,14.92_rp,-21.98_rp,   & 
   & -1.69_rp,10.68_rp,0.17_rp,0.00_rp,0.00_rp,-43.88_rp,0.00_rp,-18.93_rp,0.00_rp,   & 
   & -23.25_rp,-37.64_rp,0.00_rp,-28.65_rp,0.00_rp,0.00_rp,-10.78_rp,-23.67_rp,0.00_rp/)
    real(kind=rp), dimension(n_rows,n_cols), parameter :: table_vi_7 = (/ &     
   & 6.13_rp,8.21_rp,-27.22_rp,0.00_rp,0.00_rp,301.38_rp,0.00_rp,0.00_rp,-4.75_rp,   & 
   & -2.31_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,-59.56_rp,   & 
   & 12.66_rp,183.17_rp,-144.23_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
   & 0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
   & 3.28_rp,-15.04_rp,-3.27_rp,-56.74_rp,14.15_rp,64.74_rp,29.51_rp,201.87_rp,8.79_rp,   & 
   & 1.91_rp,17.85_rp,-45.07_rp,12.86_rp,-38.36_rp,-30.24_rp,6.90_rp,6.71_rp,-2.31_rp,   & 
   & -29.84_rp,57.93_rp,2.79_rp,0.00_rp,27.43_rp,0.00_rp,0.00_rp,0.00_rp,0.00_rp,   & 
   & 0.00_rp,-37.94_rp,102.07_rp,-13.03_rp,-134.85_rp,-69.43_rp,-21.05_rp,158.03_rp,61.37_rp,   & 
   & -1.66_rp,-11.38_rp,-4.87_rp,-20.87_rp,63.07_rp,-37.06_rp,268.47_rp,-150.70_rp,7.06_rp,   & 
   & 2.04_rp,16.31_rp,10.33_rp,-16.91_rp,2.67_rp,-18.53_rp,-5.87_rp,16.54_rp,-37.03_rp,   & 
   & 1.55_rp,37.41_rp,11.86_rp,56.87_rp,35.98_rp,-9.77_rp,-5.69_rp,17.37_rp,-60.67_rp,   & 
   & 81.75_rp,-26.52_rp,32.66_rp,-37.87_rp,-151.19_rp,29.27_rp,-15.31_rp,-1.66_rp,-73.92_rp,   & 
   & -9.46_rp,36.90_rp,-7.76_rp,165.10_rp,91.22_rp,-141.02_rp,436.49_rp,-887.25_rp,2.43_rp,   & 
   & -10.21_rp,-23.56_rp,12.97_rp,66.40_rp,28.36_rp,15.12_rp,-12.50_rp,1.82_rp,29.44_rp,   & 
   & -31.55_rp,-22.43_rp,11.33_rp,-89.17_rp,-54.51_rp,28.36_rp,-8.15_rp,-11.99_rp,60.97_rp,   & 
   & -141.60_rp,105.37_rp,-131.42_rp,33.14_rp,137.04_rp,-42.15_rp,47.86_rp,46.07_rp,40.43_rp,   & 
   & -13.73_rp,38.22_rp,-9.06_rp,229.84_rp,54.84_rp,-142.23_rp,335.42_rp,-759.51_rp,36.60_rp,   & 
   & -29.55_rp,-97.40_rp,31.14_rp,111.36_rp,129.73_rp,-77.39_rp,-12.68_rp,88.47_rp,49.66_rp,   & 
   & -12.22_rp,-85.53_rp,-13.83_rp,-95.05_rp,184.80_rp,-145.77_rp,91.04_rp,-112.43_rp,-33.31_rp,   & 
   & -218.76_rp,322.06_rp,-371.95_rp,14.40_rp,394.27_rp,50.64_rp,35.66_rp,48.81_rp,56.84_rp,   & 
   & -17.83_rp,25.34_rp,2.92_rp,206.89_rp,157.94_rp,-83.90_rp,825.13_rp,-510.85_rp,1.49_rp,   & 
   & -4.60_rp,-114.67_rp,18.73_rp,28.15_rp,-1.39_rp,-25.45_rp,-29.86_rp,-6.30_rp,65.14_rp,   & 
   & 5.36_rp,-53.26_rp,-28.29_rp,-228.34_rp,48.70_rp,-222.84_rp,226.45_rp,-59.14_rp,5.14_rp,   & 
   & -11.09_rp,28.70_rp,-28.07_rp,-20.46_rp,213.23_rp,80.24_rp,-22.57_rp,39.88_rp,84.53_rp,   & 
   & -5.16_rp,17.61_rp,47.67_rp,0.00_rp,132.35_rp,0.00_rp,0.00_rp,0.00_rp,-29.70_rp,   & 
   & 3.56_rp,-96.11_rp,0.00_rp,-13.06_rp,0.00_rp,0.00_rp,11.73_rp,-20.05_rp,-31.08_rp,   & 
   & 4.64_rp,5.13_rp,-13.92_rp,0.00_rp,0.00_rp,-122.21_rp,0.00_rp,21.94_rp,0.00_rp,   & 
   & 45.86_rp,-81.44_rp,0.00_rp,-77.55_rp,0.00_rp,0.00_rp,-50.13_rp,-78.17_rp,0.00_rp/)
    
   real(kind=rp), dimension(n_rows,n_cols,7), parameter :: table_vi(1:n_rows,1:n_cols,1:7)  = & 
      & (/table_vi_1,table_vi_2,table_vi_3,table_vi_4,table_vi_5,table_vi_6,table_vi_7/) 
    
end module choice_data
