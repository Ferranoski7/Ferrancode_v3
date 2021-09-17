import numpy as np

def rotation_z(position,angle):
    rotated=np.array([[position[i,0]*np.cos(angle)-position[i,1]*np.sin(angle),
                      position[i,0]*np.sin(angle)+position[i,1]*np.cos(angle),
                      position[i,2]] for i in range(position.shape[0])])
    return rotated

def rotation(position,axis,angle):
    c=np.cos(angle)
    s=np.sin(angle)
    value=np.sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
    x=axis[0]/value
    y=axis[1]/value
    z=axis[2]/value

    rotated=np.zeros((position.shape[0],position.shape[1]))

    for i in range(position.shape[0]):
        rotated[i,0]=(c+x*x*(1-c))*position[i,0]+(x*y*(1-c)-z*s)*position[i,1]+(x*z*(1-c)+y*s)*position[i,2]
        rotated[i,1]=(y*x*(1-c)+z*s)*position[i,0]+(c+y*y*(1-c))*position[i,1]+(y*z*(1-c)-x*s)*position[i,2]
        rotated[i,2]=(z*x*(1-c)-y*s)*position[i,0]+(z*y*(1-c)+x*s)*position[i,1]+(c+z*z*(1-c))*position[i,2]
    return rotated

