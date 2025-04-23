# in class exercise
def is_odd(i):
    '''
    Input: i, a poditive int
    Returns True if i is odd, otherwise False
    '''
    print("inside is odd")
    return i%2 != 0

print(is_odd(4))


class Coordinate (object):
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def distance(self, other):
        x_diff_sq = (self.x-other.x)


c = Coordinate(3,4)
origin = Coordinate(0,0)
print(c.x)
print(origin.y)
