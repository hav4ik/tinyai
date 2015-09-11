#include "level.h"

#include <iostream>
#include "TINY/tinyxml2.h"

namespace pl {

int Object::GetPropertyInt(std::string name){
    return atoi(properties[name].c_str());
}

float Object::GetPropertyFloat(std::string name){
    return strtod(properties[name].c_str(), NULL);
}

std::string Object::GetPropertyString(std::string name){
    return properties[name];
}

bool Level::LoadFromFile(std::string filename){
std::cout << "hello\n";

    tinyxml2::XMLDocument levelFile;

    // Загружаем XML-карту            
    tinyxml2::XMLError xmlerror = levelFile.LoadFile(filename.c_str());

    if (xmlerror == tinyxml2::XML_ERROR_FILE_NOT_FOUND){
        std::cout << "Loading level \"" << filename << "\" failed." << std::endl;
        return false;
    }
std::cout << "xml file loaded\n";

    // Работаем с контейнером map
    tinyxml2::XMLElement *map;
    map = levelFile.FirstChildElement("map");
std::cout << "map child element found\n";
    // Пример карты: <map version="1.0" orientation="orthogonal"
    // width="10" height="10" tilewidth="34" tileheight="34">
    width = atoi(map->Attribute("width"));
    height = atoi(map->Attribute("height"));
    tileWidth = atoi(map->Attribute("tilewidth"));
    tileHeight = atoi(map->Attribute("tileheight"));
std::cout << "map attribute readed" <<" "<< width<<" " << height<<" " << tileWidth<<" " << tileHeight<<" " << std::endl;
    // Берем описание тайлсета и идентификатор первого тайла
    tinyxml2::XMLElement *tilesetElement;
    tilesetElement = map->FirstChildElement("tileset");
    firstTileID = atoi(tilesetElement->Attribute("firstgid"));
std::cout << "tileset found\n";
    // source - путь до картинки в контейнере image
    tinyxml2::XMLElement *image;
    image = tilesetElement->FirstChildElement("image");
std::cout << "child element 'image' found\n";
    std::string imagepath = image->Attribute("source");
std::cout << "image source found\n";
    // Пытаемся загрузить тайлсет
    sf::Image img;
    if(!img.loadFromFile("maps/"+imagepath)){
        std::cout << "Failed to load tile sheet." << std::endl;
        return false;
    }

    // Очищаем карту от света (109, 159, 185)
    // Вообще-то в тайлсете может быть фон любого цвета, но я не нашел решения, как 16-ричную строку
    // вроде "6d9fb9" преобразовать в цвет
    img.createMaskFromColor(sf::Color(109, 159, 185));
    // Грузим текстуру из изображения
    tilesetImage.loadFromImage(img);
    // Расплывчатость запрещена
    tilesetImage.setSmooth(false);

    // Получаем количество столбцов и строк тайлсета
    int columns = tilesetImage.getSize().x / tileWidth;
    int rows = tilesetImage.getSize().y / tileHeight;

    // Вектор из прямоугольников изображений (TextureRect)
    std::vector< sf::Rect<int> > subRects;
std::cout << "creating vector\n";
    for(int y = 0; y < rows; y++){
        for(int x = 0; x < columns; x++){
            sf::Rect<int> rect;

            rect.top = y * tileHeight;
            rect.height = tileHeight;
            rect.left = x * tileWidth;
            rect.width = tileWidth;

            subRects.push_back(rect);
        }
    }


    // Работа со слоями
    tinyxml2::XMLElement *layerElement;
    layerElement = map->FirstChildElement("layer");
    while(layerElement){
        Layer layer;

        // Если присутствует opacity, то задаем прозрачность слоя, иначе он полностью непрозрачен
        if (layerElement->Attribute("opacity") != NULL){
            float opacity = strtod(layerElement->Attribute("opacity"), NULL);
            layer.opacity = 255 * opacity;
        }
        else {
            layer.opacity = 255;
        }

        // Контейнер <data>
        tinyxml2::XMLElement *layerDataElement;
        layerDataElement = layerElement->FirstChildElement("data");

        if(layerDataElement == NULL){
            std::cout << "Bad map. No layer information found." << std::endl;
        }

        // Контейнер <tile> - описание тайлов каждого слоя
        tinyxml2::XMLElement *tileElement;
        tileElement = layerDataElement->FirstChildElement("tile");

        if(tileElement == NULL){
            std::cout << "Bad map. No tile information found." << std::endl;
            return false;
        }

        int x = 0;
        int y = 0;

        while(tileElement)
        {            
            int tileGID = atoi(tileElement->Attribute("gid"));
            int subRectToUse = tileGID - firstTileID;

            // Устанавливаем TextureRect каждого тайла
            if (subRectToUse >= 0){
                sf::Sprite sprite;
                sprite.setTexture(tilesetImage);
                sprite.setTextureRect(subRects[subRectToUse]);
                sprite.setPosition(x * tileWidth, y * tileHeight);
                sprite.setColor(sf::Color(255, 255, 255, layer.opacity));

                layer.tiles.push_back(sprite);
            }

            tileElement = tileElement->NextSiblingElement("tile");

            x++;
            if (x >= width){
                x = 0;
                y++;
                if(y >= height)
                y = 0;
            }
        }

        layers.push_back(layer);

        layerElement = layerElement->NextSiblingElement("layer");        
    }

    std::cout << "level.cpp: Layer loaded\n";

    // Работа с объектами
    tinyxml2::XMLElement *objectGroupElement;

    // Если есть слои объектов
    if (map->FirstChildElement("objectgroup") != NULL)
    {
        objectGroupElement = map->FirstChildElement("objectgroup");
        while (objectGroupElement)
        {
            // Контейнер <object>
            tinyxml2::XMLElement *objectElement;
            objectElement = objectGroupElement->FirstChildElement("object");

            while(objectElement != NULL)
            {
                // Reading information about new object
                Object object;
                // Получаем все данные - тип, имя, позиция, etc
                std::string objectType;
                if (objectElement->Attribute("type") != NULL)
                {
                    objectType = objectElement->Attribute("type");
                }
                std::string objectName;
                if (objectElement->Attribute("name") != NULL)
                {
                    objectName = objectElement->Attribute("name");
                }
                int x = atoi(objectElement->Attribute("x"));
                int y = atoi(objectElement->Attribute("y"));

                if (objectType == "rect")
                {
                    int width, height;

                    sf::Sprite sprite;
                    sprite.setTexture(tilesetImage);
                    sprite.setTextureRect(sf::Rect<int>(0,0,0,0));
                    sprite.setPosition(x, y);

                    if (objectElement->Attribute("width") != NULL)
                    {
                        width = atoi(objectElement->Attribute("width"));
                        height = atoi(objectElement->Attribute("height"));
                    }
                    else
                    {
                        width = subRects[atoi(objectElement->Attribute("gid")) - firstTileID].width;
                        height = subRects[atoi(objectElement->Attribute("gid")) - firstTileID].height;
                        sprite.setTextureRect(subRects[atoi(objectElement->Attribute("gid")) - firstTileID]);
                    }

                    // Экземпляр объекта
                    object.name = objectName;
                    object.type = objectType;
                    object.sprite = sprite;

                    sf::Rect <int> objectRect;
                    objectRect.top = y;
                    objectRect.left = x;
                    objectRect.height = height;
                    objectRect.width = width;
                    object.rect = objectRect;

                    // "Переменные" объекта
                    tinyxml2::XMLElement *properties;
                    properties = objectElement->FirstChildElement("properties");
                    if (properties != NULL)
                    {
                        tinyxml2::XMLElement *prop;
                        prop = properties->FirstChildElement("property");
                        if (prop != NULL)
                        {
                            while(prop != NULL)
                            {
                                std::string propertyName = prop->Attribute("name");
                                std::string propertyValue = prop->Attribute("value");

                                object.properties[propertyName] = propertyValue;

                                prop = prop->NextSiblingElement("property");
                            }
                        }
                    }
                }
                if (objectType == "polygon")
                {
                    // Экземпляр объекта
                    object.name = objectName;
                    object.type = objectType;

                    // points of this polygon
                    tinyxml2::XMLElement *points;
                    std::string pointArray;
                    points = objectElement->FirstChildElement("polygon");
                    if (points != NULL)
                    {
                        pointArray = points->Attribute("points");
                    }
                    std::cout << "points: " << pointArray << std::endl;
                    for (int h=0; h<pointArray.length(); h++)
                    if (pointArray[h] == ',')
                    pointArray[h] = ' ';
                    std::cout << pointArray << std::endl;
                    std::stringstream ss(pointArray);

                    do
                    {
                        // read as many numbers as possible.
                        for (int tx, ty; ss >> tx >> ty;) {
                            object.points.push_back(sf::Vector2i(x+tx,y+ty));
                        }
                        // consume and discard token from stream.
                        if (ss.fail())
                        {
                            ss.clear();
                            std::string token;
                            ss >> token;
                        }
                    }
                    while (!ss.eof());

                    // "Переменные" объекта
                    tinyxml2::XMLElement *properties;
                    properties = objectElement->FirstChildElement("properties");
                    if (properties != NULL)
                    {
                        tinyxml2::XMLElement *prop;
                        prop = properties->FirstChildElement("property");
                        if (prop != NULL)
                        {
                            while(prop != NULL)
                            {
                                std::string propertyName = prop->Attribute("name");
                                std::string propertyValue = prop->Attribute("value");

                                object.properties[propertyName] = propertyValue;

                                prop = prop->NextSiblingElement("property");
                            }
                        }
                    }
                }
                // Пихаем объект в вектор
                objects.push_back(object);

                objectElement = objectElement->NextSiblingElement("object");
            }
            objectGroupElement = objectGroupElement->NextSiblingElement("objectgroup");
        }
    }
    else
    {
        std::cout << "No object layers found..." << std::endl;
    }

    return true;
}

Object Level::GetObject(std::string name){
    // Только первый объект с заданным именем

    for (int i = 0; i < objects.size(); i++)
        if (objects[i].name == name)
            return objects[i];
}

std::vector<Object> Level::GetObjects(std::string name)
{
    // Все объекты с заданным именем
    std::vector<Object> vec;

    for(int i = 0; i < objects.size(); i++)
        if(objects[i].name == name)
            vec.push_back(objects[i]);

    return vec;
}

sf::Vector2i Level::GetTileSize()
{
    return sf::Vector2i(tileWidth, tileHeight);
}

void Level::Draw(sf::RenderWindow &window)
{
    // Рисуем все тайлы (объекты НЕ рисуем!)
    for(int layer = 0; layer < layers.size(); layer++)
        for(int tile = 0; tile < layers[layer].tiles.size(); tile++)
            window.draw(layers[layer].tiles[tile]);
}
int Level::GetHeight(){
    return height * tileHeight;
}
int Level::GetWidth(){
    return width * tileWidth;
}

} // end of namespace pl
