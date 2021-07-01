

class BoundaryCell{
public:
    BoundaryCell(float x, float y, bool left, bool top, bool right, bool bottom);
    void draw();

    float b_x;
    float b_y;
    bool b_left;
    bool b_top;
    bool b_right;
    bool b_bottom;
};