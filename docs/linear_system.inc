.. math::

    \left[ \begin{array}{ccccc}
    A^0_1      & A^{1}_1    & \textbf{0} & \cdots       & \textbf{0}    & A^{-1}_1   \\
    A^{-1}_2   & A^0_2      & \ddots     & \ddots       &               & \textbf{0} \\
    \textbf{0} & \ddots     & \ddots     & \ddots       & \ddots        & \vdots     \\
    \vdots     & \ddots     & \ddots     & \ddots       & \ddots        & \textbf{0} \\
    \textbf{0} &            & \ddots     & \ddots       & A^0_{ny}      & A^1_{ny}   \\
    A^1_{ny+1} & \textbf{0} & \cdots     & \textbf{0}   & A^{-1}_{ny+1} & A^0_{ny+1} \\
    \end{array} \right]
    \left[ \begin{array}{c}
    x_1 \\ x_2 \\ \vdots \\ \vdots \\ x_{ny} \\ x_{ny+1} \\
    \end{array} \right]
    =
    \left[ \begin{array}{c}
    y_1 \\ y_2 \\ \vdots \\ \vdots \\ y_{ny} \\ y_{ny+1} \\
    \end{array} \right]
    
    \text{where}
    
    A^k_i =
    \left\{ \begin{array}{rl}
    \left[ \begin{array}{ccccccc}
    1              & -2               & 1                & 0
                   & \cdots           & \cdots           & 0                 \\
    a^{0,-1}_{i,2} & a^{0,0}_{i,2}    & a^{0,1}_{i,2}    & \ddots
                   &                  &                  & \vdots            \\
    0              & a^{0,-1}_{i,3}   & a^{0,0}_{i,3}    & a^{0,1}_{i,3}
                   & \ddots           &                  & \vdots            \\
    \vdots         & \ddots           & \ddots           & \ddots
                   & \ddots           & \ddots           & \vdots            \\
    \vdots         &                  & \ddots           & a^{0,-1}_{i,nx-1}
                   & a^{0,0}_{i,nx-1} & a^{0,1}_{i,nx-1} & 0                 \\
    \vdots         &                  &                  & 0
                   & a^{0,-1}_{i,nx}  & a^{0,0}_{i,nx}   & a^{0,1}_{i,nx}    \\
    0              & \cdots           & \cdots           & 0
                   & 1                & -2               & 1                 \\
    \end{array} \right] & \text{if } k = 0 \\
    \left[ \begin{array}{ccccccc}
    0               & 0                & \cdots           & \cdots
                    & \cdots           & \cdots           & 0                 \\
    0               & a^{k,0}_{i,2}    & \ddots           &
                    &                  &                  & \vdots            \\
    \vdots          & \ddots           & a^{k,0}_{i,3}    & \ddots
                    &                  &                  & \vdots            \\
    \vdots          &                  & \ddots           & \ddots
                    & \ddots           &                  & \vdots            \\
    \vdots          &                  &                  & \ddots
                    & a^{k,0}_{i,nx-1} & \ddots           & \vdots            \\
    \vdots          &                  &                  &
                    & \ddots           & a^{k,0}_{i,nx}   & 0                 \\
    0               & \cdots           & \cdots           & \cdots  
                    & \cdots           & 0                & 0                 \\
    \end{array} \right] & \text{if } k \neq 0 \\
    \end{array} \right.
  
    \begin{array}{rclrcl}
    x_i &=& 
    \left[ \begin{array}{c}
    C_{t,i,1}^{n+1} \\
    C_{t,i,2}^{n+1} \\
    C_{t,i,3}^{n+1} \\
    \vdots \\
    C_{t,i,nx-1}^{n+1} \\
    C_{t,i,nx}^{n+1} \\
    C_{t,i,nx+1}^{n+1} \\
    \end{array} \right] &
    y_i &=&
    \left[ \begin{array}{c}
    0 \\
    y_{i,2} \\
    y_{i,3} \\
    \vdots \\
    y_{i,nx-1} \\
    y_{i,nx} \\
    0 \\
    \end{array} \right] \\
    \end{array}
     
    \text{and where}
     
    \begin{array}{rcll}
    a^{0,0}_{i,j} &=& 1 + \alpha \left( \delta_x u_{w,x,i,j} \frac{\Delta t}{\Delta x_{i,j}} + \delta_y u_{w,y,i,j} \frac{\Delta t}{\Delta y_{i,j}} + \frac{\Delta t}{T} \right) \\
    a^{-1,0}_{i,j} &=& -\alpha u_{w,y,i,j} \frac{\Delta t}{\Delta y_i,j} \beta \\
    a^{0,-1}_{i,j} &=& -\alpha u_{w,x,i,j} \frac{\Delta t}{\Delta x_i,j} \beta \\
    a^{0,1}_{i,j} &=& \alpha u_{w,x,i,j} \frac{\Delta t}{\Delta x_i,j} \left( 1 - \beta \right) \\
    a^{1,0}_{i,j} &=& \alpha u_{w,y,i,j} \frac{\Delta t}{\Delta y_i,j} \left( 1 - \beta \right) \\
    y_{i,j} &=& C^n_{t,i,j} + \alpha w_{i,j} \frac{\Delta t}{T} C_{u,i,j} + (1 - \alpha) \left[ \phantom{\frac{\Delta t}{T}} \right. & w_{0,i,j} \frac{\Delta t}{T} C_{u,i,j} \\
      & & & - \left( \delta_x u_{w,x,i,j} \frac{\Delta t}{\Delta x_{i,j}} + \delta_y u_{w,y,i,j} \frac{\Delta t}{\Delta y_{i,j}} + \frac{\Delta t}{T} \right) C_{t,i,j} \\
      & & & + u_{w,x,i,j} \frac{\Delta t}{\Delta x_{i,j}} \beta C_{t,i-1,j} \\
      & & & - u_{w,x,i,j} \frac{\Delta t}{\Delta x_{i,j}} \left( 1 - \beta \right) C_{t,i+1,j} \\
      & & & + u_{w,y,i,j} \frac{\Delta t}{\Delta y_{i,j}} \beta C_{t,i,j-1} \\
      & & & - u_{w,y,i,j} \frac{\Delta t}{\Delta y_{i,j}} \left( 1 - \beta \right) C_{t,i,j+1} \left. \phantom{\frac{\Delta t}{T}} \right] \\
    \end{array}

