class BaseObject
{
public:
    BaseObject();

    virtual void draw() = 0;
    virtual void update() = 0;
    virtual void setPosition(float x, float y) = 0;
    virtual void setVelocity(float u, float v) = 0;
    virtual bool isOnCell(float x, float y) = 0;
};