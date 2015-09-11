#include <iostream>
#include <fstream>

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>

#include <vector>
#include <cmath>

#include "TINY/tinyxml2.h"
#include "level.h"
#include "physics.h"

#include <tinyann.hpp>
#include <tinyneat.hpp>

using namespace std;

class car {
    sf::Vector2f spawn_point;
    sf::Rect<int> destination;

    sf::Vector2f direction;
    sf::Vector2f speed;
    sf::Vector2f position;

    sf::Vector2f sensor_c;
    sf::Vector2f sensor_l;
    sf::Vector2f sensor_ll;
    sf::Vector2f sensor_r;
    sf::Vector2f sensor_rr;

    sf::Clock cl;
    float rotation; // positive - left, negative - right, gradus per sec

    unsigned int fitness;
    float stale_time;
    bool dead;
public:
    unsigned int get_fitness(){
        return fitness;
    }

    bool is_winner(){
        bool inside = true;
        inside = inside && (position.x >= destination.left) && (position.x <= destination.left + destination.width);
        inside = inside && (position.y >= destination.top) && (position.y <= destination.top + destination.height);
        if (inside){
            float w = spawn_point.x - destination.left;
            float h = spawn_point.y - destination.top;
            fitness = std::sqrt(w*w + h*h) * 100;
        }
        return inside;
    }

    car (sf::Vector2f pos, sf::Rect<int> finish) {
        dead = false;
        fitness = 0;
        stale_time = 0.0;

        rotation = 0.0; // gradus per sec
        spawn_point = pos;
        destination = finish;
        position = pos;
        direction = sf::Vector2f(0.0, -1.0);

        sensor_c.x = direction.x * 50.0;
        sensor_c.y = direction.y * 50.0;

        float rl = -30.0 / 180.0 * M_PI;
        float rll = -60.0 / 180.0 * M_PI;
        float rr = 30.0 / 180.0 * M_PI;
        float rrr = 60.0 / 180.0 * M_PI;

        sensor_l  = sf::Vector2f(sensor_c.x * cos(rl) - sensor_c.y * sin(rl), sensor_c.x * sin(rl) + sensor_c.y * cos(rl));
        sensor_ll = sf::Vector2f(sensor_c.x * cos(rll) - sensor_c.y * sin(rll), sensor_c.x * sin(rll) + sensor_c.y * cos(rll));
        sensor_r  = sf::Vector2f(sensor_c.x * cos(rr) - sensor_c.y * sin(rr), sensor_c.x * sin(rr) + sensor_c.y * cos(rr));
        sensor_rr = sf::Vector2f(sensor_c.x * cos(rrr) - sensor_c.y * sin(rrr), sensor_c.x * sin(rrr) + sensor_c.y * cos(rrr));

    }

    void rotate(float angle){
        float rad = M_PI * angle / 180.0;
        float x, y;
        x = direction.x * cos(rad) - direction.y * sin(rad);
        y = direction.x * sin(rad) + direction.y * cos(rad);
        direction = sf::Vector2f(x, y);

        sensor_c.x = direction.x * 50.0;
        sensor_c.y = direction.y * 50.0;

        float rl = -30.0 / 180.0 * M_PI;
        float rll = -60.0 / 180.0 * M_PI;
        float rr = 30.0 / 180.0 * M_PI;
        float rrr = 60.0 / 180.0 * M_PI;

        sensor_l  = sf::Vector2f(sensor_c.x * cos(rl) - sensor_c.y * sin(rl), sensor_c.x * sin(rl) + sensor_c.y * cos(rl));
        sensor_ll = sf::Vector2f(sensor_c.x * cos(rll) - sensor_c.y * sin(rll), sensor_c.x * sin(rll) + sensor_c.y * cos(rll));
        sensor_r  = sf::Vector2f(sensor_c.x * cos(rr) - sensor_c.y * sin(rr), sensor_c.x * sin(rr) + sensor_c.y * cos(rr));
        sensor_rr = sf::Vector2f(sensor_c.x * cos(rrr) - sensor_c.y * sin(rrr), sensor_c.x * sin(rrr) + sensor_c.y * cos(rrr));

        rotation = 0.0;
    }

    void move(std::vector<ph::intPolygon>& pol){
        if (dead)
            return ;

        float elapsed = cl.restart().asSeconds();

        float angle = rotation * elapsed;
        rotate(angle);

        sf::Vector2f v;
        v.x = speed.x * elapsed;
        v.y = speed.y * elapsed;

        sf::Vector2f np;
        np.x = position.x + v.x;
        np.y = position.y + v.y;

        for (size_t i=0; i<pol.size(); i++){
            ph::intPolygon& p = pol[i];
            sf::Vector2f ip = p.CheckIntersect(position, position + v);
            if (ip.x >= 0 && ip.y >= 0)
                np = sf::Vector2f(ip.x - v.x/100.0, ip.y - v.y/100.0);
        }

        speed.x = 0.0;
        speed.y = 0.0;
        position = np;

        float a = position.x - spawn_point.x;
        float b = position.y - spawn_point.y;
        unsigned int new_fitness = std::sqrt(a*a + b*b);
        if (new_fitness > fitness){
            stale_time = 0.0;
            fitness = new_fitness;
        }

        stale_time += elapsed;
        if (stale_time > 5.0)
            dead = true;

    }

    void rotate_right(){
        rotation = 100.0;
    }

    void rotate_left(){
        rotation = -100.0;
    }

    void move_forward(){
        speed.x = direction.x * 100;
        speed.y = direction.y * 100;
    }

    void Draw(sf::RenderWindow& window){
        sf::CircleShape c;
        c.setFillColor(sf::Color::Black);
        c.setRadius(8.0);
        c.setPosition(position.x - 8, position.y - 8);

        sf::Vertex line1[] = { sf::Vertex(position), sf::Vertex(position + sensor_c) };
        window.draw(line1, 2, sf::Lines);

        sf::Vertex line2[] = { sf::Vertex(position), sf::Vertex(position + sensor_l) };
        window.draw(line2, 2, sf::Lines);

        sf::Vertex line3[] = { sf::Vertex(position), sf::Vertex(position + sensor_ll) };
        window.draw(line3, 2, sf::Lines);

        sf::Vertex line4[] = { sf::Vertex(position), sf::Vertex(position + sensor_r) };
        window.draw(line4, 2, sf::Lines);

        sf::Vertex line5[] = { sf::Vertex(position), sf::Vertex(position + sensor_rr) };
        window.draw(line5, 2, sf::Lines);



        window.draw(c);
    }
};

void output_info(
        sf::RenderWindow& window,
        sf::Font& font,
        unsigned int generation,
        unsigned int specie,
        unsigned int genomes,
        unsigned int global_maxfitness,
        unsigned int current_specie_max_fitness,
        bool have_winner
    ){

    sf::Color bgcolor;
    bgcolor.g = 0;
    bgcolor.r = 0;
    bgcolor.b = 0;
    bgcolor.a = 162;

    sf::RectangleShape background;
    background.setPosition(window.getSize().x - 48*6, window.getSize().y - 48*3);
    background.setSize(sf::Vector2f(48*6, 48*3));
    background.setFillColor(bgcolor);
    background.setOutlineThickness(3);
    background.setOutlineColor(sf::Color::Black);

    std::stringstream ss;
    ss << "Generation: " << generation << std::endl;
    ss << "Specie number: " << specie << std::endl;
    ss << "Genomes in specie: " << genomes << std::endl;
    ss << "Global max fitness: " << global_maxfitness << std::endl;
    ss << "Current specie max: " << current_specie_max_fitness << std::endl;
    if (have_winner)
        ss << "WE HAVE A WINNER!!!" << std::endl;

    sf::Text text;
    text.setFont(font);
    text.setPosition(sf::Vector2f(background.getPosition().x+6, background.getPosition().y+2));
    text.setString(ss.str());
    text.setCharacterSize(20);
    text.setColor(sf::Color::Yellow);

    window.draw(background);
    window.draw(text);
}

int main()
{
    pl::Level level;
    level.LoadFromFile("maps/track.tmx");
    sf::Font font;
    font.loadFromFile("UbuntuMono-R.ttf");

    sf::RenderWindow window;
    window.create(sf::VideoMode(768, 768), "Kolobosha adventures");
    window.setVerticalSyncEnabled(true); // call it once, after creating the window
    window.setFramerateLimit(30.0); // call it once, after creating the window

    pl::Object start = level.GetObject("start");
    pl::Object finish = level.GetObject("finish");
    std::vector<pl::Object> obj = level.GetObjects("wall");    

    std::vector<ph::intPolygon> pol;
    for (size_t i=0; i<obj.size(); i++){
        ph::intPolygon p(obj[i].rect);
        pol.push_back(p);
    }



    car c(sf::Vector2f(start.rect.left, start.rect.top), finish.rect);

    

    while(window.isOpen())
    {
        sf::Keyboard key;
        if (key.isKeyPressed(sf::Keyboard::A))
            c.rotate_left();
        if (key.isKeyPressed(sf::Keyboard::D))
            c.rotate_right();
        if (key.isKeyPressed(sf::Keyboard::W))
            c.move_forward();

        sf::Event event;
        while(window.pollEvent(event))
        {
            if(event.type == sf::Event::Closed)
                window.close();

        }

        window.clear();

        level.Draw(window);
        c.move(pol);
        c.Draw(window);
        output_info(window, font, 1, 1, 10, c.get_fitness(), c.get_fitness(), c.is_winner());

        window.display();
    }

    return 0;
}

