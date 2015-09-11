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

    bool is_alive(){
        return !dead;
    }

    void get_sensor(std::vector<double>& v, std::vector<ph::intPolygon>& pol){
        sf::Vector2f sc, sl, sll, sr, srr;
        float c, l, ll, r, rr;

        sc = position + sensor_c;
        sl = position + sensor_l;
        sll = position + sensor_ll;
        sr = position + sensor_r;
        srr = position + sensor_rr;

        for (size_t i=0; i<pol.size(); i++){
            ph::intPolygon& p = pol[i];
            sf::Vector2f ip;

            ip = p.CheckIntersect(position, position + sensor_c);
            if (ip.x >= 0 && ip.y >= 0){
                float a = ip.x - position.x;
                float b = ip.y - position.y;
                c = std::sqrt(a*a + b*b);
            }

            ip = p.CheckIntersect(position, position + sensor_l);
            if (ip.x >= 0 && ip.y >= 0){
                float a = ip.x - position.x;
                float b = ip.y - position.y;
                l = std::sqrt(a*a + b*b);
            }

            ip = p.CheckIntersect(position, position + sensor_ll);
            if (ip.x >= 0 && ip.y >= 0){
                float a = ip.x - position.x;
                float b = ip.y - position.y;
                ll = std::sqrt(a*a + b*b);
            }

            ip = p.CheckIntersect(position, position + sensor_r);
            if (ip.x >= 0 && ip.y >= 0){
                float a = ip.x - position.x;
                float b = ip.y - position.y;
                r = std::sqrt(a*a + b*b);
            }

            ip = p.CheckIntersect(position, position + sensor_rr);
            if (ip.x >= 0 && ip.y >= 0){
                float a = ip.x - position.x;
                float b = ip.y - position.y;
                rr = std::sqrt(a*a + b*b);
            }                
        }

        c = c / 25.0;
        l = l / 25.0;
        ll = ll / 25.0;
        r = r / 25.0;
        rr = rr / 25.0;

        v[0] = c; v[1] = l; v[2] = ll; v[3] = r; v[4] = rr;
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
        rotation = 110.0;
    }

    void rotate_left(){
        rotation = -110.0;
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

    
    // 5 input, 3 output, 1 bias, can be recurrent
    neat::pool p(5, 3, 1, true);
    p.import_fromfile("generation.dat");
    bool have_a_winner = false;
    unsigned int global_maxfitness = 0;

    // vector of cars
    std::vector< std::pair<car,ann::neuralnet> > cars;

    // iterator
    unsigned int specie_counter = 0;
    auto specie_it = p.species.begin();

    // init initial
    if (specie_it != p.species.end())
        for (size_t i=0; i<(*specie_it).genomes.size(); i++){
            car new_car(sf::Vector2f(start.rect.left, start.rect.top), 
                finish.rect);
            ann::neuralnet n;
            n.from_genome((*specie_it).genomes[i]);
            cars.push_back(std::make_pair(new_car, n));
        }            

    // main loop
    while(window.isOpen() && (!have_a_winner))
    {        
        sf::Event event;
        while(window.pollEvent(event))
        {
            if(event.type == sf::Event::Closed)
                window.close();
            else 
                if (event.type == sf::Event::KeyPressed)
                    if (event.key.code == sf::Keyboard::Escape)
                        return 0;

        }


        bool all_dead = true;
        for (size_t i=0; i<cars.size(); i++)
            if (cars[i].first.is_alive())
                all_dead = false;

        if (all_dead){
            
            if (specie_it != p.species.end()){
                int best_id = -1;
                for (size_t i=0; i<(*specie_it).genomes.size(); i++){
                    (*specie_it).genomes[i].fitness = 
                            cars[i].first.get_fitness();
                    if ((*specie_it).genomes[i].fitness > global_maxfitness){
                    global_maxfitness = (*specie_it).genomes[i].fitness; 
                        best_id = i;
                    }       
                }
                if (best_id != -1){
                    ann::neuralnet& n = cars[best_id].second;
                    n.export_tofile("best_network");
                }
            }

            specie_it++;
            specie_counter++;
            cars.clear();

            if (specie_it == p.species.end()){                
                p.new_generation();
                std::string fname = "result/gen";
                fname += std::to_string(p.generation());
                p.export_tofile(fname);
                p.export_tofile("generation.dat");
                std::cerr << "Starting new generation. Number = " << p.generation() << std::endl;
                specie_it = p.species.begin();
                specie_counter = 0;
            }

            if (specie_it != p.species.end())
                for (size_t i=0; i<(*specie_it).genomes.size(); i++){
                    car new_car(sf::Vector2f(start.rect.left, start.rect.top), 
                        finish.rect);
                    ann::neuralnet n;
                    n.from_genome((*specie_it).genomes[i]);
                    cars.push_back(std::make_pair(new_car, n));
                }            
        }

        for (size_t i=0; i<cars.size(); i++){            
            std::vector<double> input(5, 0.0);
            std::vector<double> output(3, 0.0);

            car& c = cars[i].first;
            ann::neuralnet& n = cars[i].second;

            if (!c.is_alive())
                continue;

            c.get_sensor(input, pol);
            n.evaluate(input, output);

            if (output[0] > 0.0)
                c.move_forward();
            if (output[1] > 0.0)
                c.rotate_left();
            if (output[2] > 0.0)
                c.rotate_right();

            c.move(pol);
        }



        unsigned int local_maxfitness = 0;
        for (size_t i=0; i<cars.size(); i++)
            if (cars[i].first.get_fitness() > local_maxfitness)
                local_maxfitness = cars[i].first.get_fitness();
    

        size_t winner_id;
        for (size_t i=0; i<cars.size(); i++)
            if (cars[i].first.is_winner()){
                have_a_winner = true;
                winner_id = i;
            }

        if (have_a_winner)
            cars[winner_id].second.export_tofile("winner_network");


        window.clear();

        level.Draw(window);
        for (size_t i=0; i<cars.size(); i++)
            cars[i].first.Draw(window);

        output_info(window, font, p.generation(), 
            specie_counter, cars.size(), global_maxfitness, local_maxfitness,
            have_a_winner);
       
        window.display();
    }

    
    return 0;
}

