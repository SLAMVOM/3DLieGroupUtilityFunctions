# This repository includes the C++ implementation of several identities and functions on SO(3) and SE(3)



-----------------------------------------------

For a more thorough treatment on 3D Special Orthogonal Group SO(3) and Special Euclidean Group SE(3), one may refer to:
- [1] Prof. Barfoot's book: [State Estimation For Robotics](https://doi.org/10.1017/9781316671528)
- [2] Prof. Chirikjian's book: [Stochastic Models, Information Theory, and Lie Groups, Volume 1](https://doi.org/10.1007/978-0-8176-4803-9)
- [3] Sola et al., (2018) [A micro Lie theory for state estimation in robotics](
https://doi.org/10.48550/arXiv.1812.01537)

---------------------------------------------

#### To use the code, simply copy and paste the file into the code directory and include the file to the code. E.g.,
`
	#include "lie_utils.h"
`



The only dependency requried is the [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page) package, which a commonly used Linear Algebra pacakge in C++.

Other useful references:

- [Google Ceres-solver](http://ceres-solver.org/): the "rotation.h" file provides good tamplated functions for operations on SO(3)
- [Dr. Lee Clement's](https://github.com/utiasSTARS/liegroups) Python Implementation of Lie Group operations on GitHub
- [Dr. Martin Brossard's](https://github.com/mbrossar/SE2-3-) Python implementation of Lie Group utility and plotting functions
- and many others

--------------------------------------
The notation and convention style adopted in this repo mainly follow the one specified in [1]. 

| Symbol | Description |
| :------: | :---- |
| $\bf R$ | 3 x 3 rotation matrix on SO(3) |
| $\bf q$ | 4 x 1 unit quaternion, where  $\bf q =$  $[ q_w\ q_x\ q_y\ q_z]^{\it T}$, and $q_{w}$ is the scalar part of the unit quaternion |
| $\pmb{\phi}$ | 3 x 1 axis-angle representation of rotation |
| $\phi$ | scalar rotation angle as $\| \pmb{\phi} \|_2$ |
| $\bf a$ | 3 x 1 axis of rotation, where $\phi\bf a = \pmb{\phi}$ and $\bf a^{\it T}\bf a \rm = \rm 1$ |
| $\bf J$ | 3 x 3 left Jacobian of SO(3)|
| $\bf T$ | 4 x 4 transformation matrix on SE(3), where $\bf T \rm \equiv \begin{bmatrix} \bf R & \bf J\pmb{\rho} \\ \bf 0^{\it T} & 1 \end{bmatrix}= \begin{bmatrix} \bf R & \bf t \\ \bf 0^{\it T} & 1 \end{bmatrix} $ |
| $\bf t$ | 3 x 1 translational vector, and $\bf t = \bf J\pmb{\rho}$ |
| $\pmb{\xi}$ or $\pmb{\zeta}$ | 6 x 1 tangential vector corresponding to the transformation on SE(3), where $$\pmb{\xi} = \begin{bmatrix} \pmb \rho \\ \pmb \phi\end{bmatrix}$$ |
| $ \mathcal{T}$ | 6 x 6 adjoint form of transformation matrix, where  $\mathcal{T} = \rm Ad(\bf T) = \begin{bmatrix} \bf R & (\bf J \pmb{\rho})^{\wedge}\bf R \\ \bf 0_{\small 3\times3} & \bf R\end{bmatrix}$ |
| $\mathcal{J}$ | 6 x 6 left Jacobian of SE(3) |
| $\cdot^{\wedge}$ | skew-symmetric operator (overloaded), $$ \pmb{a}^{\wedge} = \begin{bmatrix} a_1 \\ a_2 \\ a_3 \end{bmatrix}^{\wedge} = \begin{bmatrix} 0 & -a_3 & a_2 \\ a_3 & 0 & -a_1 \\ -a_2 & a_1 & 0 \end{bmatrix}_{3 \times 3} $$  and  $$\pmb \xi^{\wedge} = \begin{bmatrix} \pmb \rho \\ \pmb \phi \end{bmatrix}^{\wedge} = \begin{bmatrix} \pmb \phi^{\wedge} & \pmb \rho \\ \bf 0^{\it T} & 0 \end{bmatrix}_{4 \times 4}$$ |
| $\cdot^{\curlywedge}$ | curly hat operator over SE(3), $$\pmb \xi^{\curlywedge} = \begin{bmatrix} \pmb \rho \\ \pmb \phi \end{bmatrix}^{\curlywedge} = \begin{bmatrix} \pmb \phi^{\wedge} & \pmb \rho^{\wedge} \\ \bf 0_{\rm \small 3\times3} & \pmb \phi^{\wedge} \end{bmatrix}_{6 \times 6}$$ |
| $\rm exp(\cdot)$ | matrix exponential map (overloaded), where  $\bf R = \rm exp(\pmb \phi^{\wedge})$  and  $\bf T = \rm exp(\pmb \xi^{\wedge}) $  and  $\mathcal{T} = \rm exp(\pmb \xi^{\curlywedge})$ |



SO(3) utilities:

| Function | Usage | Formula / Identity |
| ---- | ---- | ---- |
| unitQuaternionToAngleAxis(T* q, T* aa)<br /><br />This is a template function.<br />**Input**: `q` - 1D array with 4 elements of quaternion<br />**Output**: `aa` - 1D array with 3 elements of angle-axis | `double aa[3];`<br />`// input q is quaternion`<br />`unitQuaternionToAngleAxis(q, aa);` | $ k = \begin{cases} \rm{atan2} (\sqrt{ q_x^2 + q_y^2 + q_z^2}\ ,\ q_w)\ & \rm if\  q_w\lt0 \\ \rm{atan2} (-\sqrt{ q_x^2 + q_y^2 + q_z^2}\ ,-q_w)\ &\rm if\  q_w\lt0 \end{cases} \\ \\ \pmb \phi = \begin{bmatrix} \phi_{1} \\ \phi_{2} \\ \phi_{3} \end{bmatrix} = \begin{bmatrix} q_x \frac{2k}{\sqrt{ q_x^2 + q_y^2 + q_z^2}} \\ q_y \frac{2k}{\sqrt{ q_x^2 + q_y^2 + q_z^2}} \\ q_z \frac{2k}{\sqrt{ q_x^2 + q_y^2 + q_z^2}} \end{bmatrix}$ |
| angleAxisToUnitQuaternion(T* aa, T* q)<br /><br /> This is a template function.<br />**Input**: `aa` - 1D array with 3 elements of angle-axis<br />**Output**: `q` - 1D array with 4 elements of quaternion |`double q[4];`<br />`// input aa is angle-axis`<br />`angleAxisToUnitQuaternion(aa, q);` | $\phi = \|\pmb \phi\|,\quad \bf{a} = \frac{\pmb\phi}{\phi} \\ \begin{bmatrix} \ q_w \\ q_x \\ q_y \\ q_z \end{bmatrix} = \begin{bmatrix} \cos(\frac{\phi}{2}) \\ a_1\sin(\frac{\phi}{2}) \\ a_2\sin(\frac{\phi}{2}) \\ a_3\sin(\frac{\phi}{2}) \end{bmatrix}$ |
| rotationMatrixArrayToUnitQuaternion(T R_[9], T q[4])<br /><br />This is a template function.<br />**Input**: `R_` - 1D array |      |      |









------------------------------------------------

#### Author: MT
#### Creation Date: 2022-May-08
#### Previous Edit: 2022-May-08
