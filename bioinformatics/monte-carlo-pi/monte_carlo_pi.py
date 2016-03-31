from random import random

def monte_carlo_pi(intCount):
    count_circle = 0
    count_rect = 0

    for i in range(intCount):
        rand_x = random()
        rand_y = random()

        if rand_x**2 + rand_y**2 < 1:
            count_circle += 1
        count_rect += 1

    return 4 * count_circle / count_rect

print(monte_carlo_pi(1000000))