from numpy import real


class Complex:
    def __init__(self,realpart,imagpart):
        self.r=realpart
        self.i=imagpart
x=Complex(2,3)
print('{}+{}i'.format(x.r,x.i))