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

Some notes regarding the usage of the file, please refer to [3D_Utilities_Note.pdf](3D_Utilities_Note.pdf) in this repository.

<br>
<br>
<br>
<br>

Other useful references:
- [Google Ceres-solver](http://ceres-solver.org/): the "rotation.h" file provides good tamplated functions for operations on SO(3)
- [Dr. Lee Clement's](https://github.com/utiasSTARS/liegroups) Python Implementation of Lie Group operations on GitHub
- [Dr. Martin Brossard's](https://github.com/mbrossar/SE2-3-) Python implementation of Lie Group utility and plotting functions
- and many others




------------------------------------------------

#### Author: MT
#### Creation Date: 2022-May-08
#### Previous Edit: 2022-May-09
