# Example of simple phylogenetic tree drawn by Python module called `turtle`.

Tree is drawn recursively, first are left branches and then right branches.
The idea is simple: `turtle` always stay on its last position and direction.
So we start heading to the east (0°) -- this is our default direction.
Then we draw the left branches until ending node appears. When it happens, we draw
the ending node label, move back and set our direction to default (east, 0°).
Then we draw the right branches similar to left branches, except we use opposite angles.
The principle is in returning of turtle back to its default direction.

Just remind that branches are not proportional as should be!

# Example of recursive drawing of tree.

The second example is showing the recursive drawing of tree, which stops after
some conditions are met (`branch length < min` etc.). You can edit this script
to draw a tree more randomly/naturally - randomize branch angles, farther branches
are thinner, etc.