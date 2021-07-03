#include "Square.h"
#include "gfx/mat2.h"

class SquareRigid : BaseObject
{
public:
    /**
         * @param x: x coordinate of the center
         * @param y: y coordinate of the center
         * @param height: height of the object
         * @param width: width of the object
         * @param N: grid size
         */
    SquareRigid(float x, float y, int height, int width, int N);
    void draw() override;
    void update() override;
    void setPosition(float x, float y) override;
    void setVelocity(float u, float v) override;
    bool isOnCell(float x, float y) override;
    std::vector<float> getPosition();
    std::vector<float> getVelocity();
    std::vector<float> getBoundaries();
    
    Mat2 R;
    float omega;
    float torque;
    float M;
    std::vector<Vec2f> F;
    Mat2 Ibody, Ibodyinv;

private:
    float o_x;
    float o_y;
    float o_u;
    float o_v;
    int o_h;
    int o_w;
    int grid_size;
};