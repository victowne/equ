# equ
=-------------------------------------=  
 Author: Weikang Tang (Vic Towne)    
 Mail: kanger@mail.dlut.edu.cn   
=-------------------------------------=   

compatible with Python 3.7 

libraries required:  
numpy,h5py,matplotlib,scipy

reqdsk.py -> read geqdsk file format & interpolation  
poincareRK4.py -> equlibrium magnetic field & creat triangle mesh  

try:  
python main.py gfilename nlayers  

Some examples:  
![Bphi](/picture/Bphi.png)  
Fig 1. contour plot of equilibrium toroidal magnetic field  
![poincare](/picture/pc.png)  
Fig 2. poincare plot of magnetic field  
![triagnle](/picture/mesh.png)  
Fig 3. triangle mesh based on the magnetic flux  
