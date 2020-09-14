//  Student.cpp

#include "student.h"

// Constructor
Student::Student(std::string name, int age, bool male) {
    this->name = name;
    this->age = age;
    this->male = male;
    this->favoriteNumbers = {2, 3, 5, 7, 11};
}

// Getters
bool Student::IsMale() { return male; }
int Student::GetAge() { return age; }
std::string Student::GetName() { return name; }
std::vector<int> Student::GetFavoriteNumbers() { return favoriteNumbers; }
void Student::SetFavoriteNumbers(const std::vector<int>& data) { favoriteNumbers = data; }
std::vector<int> Student::GetNumbers() { return {1, 2, 3}; }

BadStudent::BadStudent(std::string name, int age, bool male) : Student(name, age, male) {}
std::vector<int> BadStudent::GetNumbers() { return {-1, -2, -3}; }

// Methods
bool Student::LikesBlue() { return (male || age >= 10); }

Teacher::Teacher(std::string name, int age, bool male) {
    this->name = name;
    this->age = age;
    this->male = male;
}

// Getters
bool Teacher::IsMale() { return male; }
int Teacher::GetAge() { return age; }
std::string Teacher::GetName() { return name; }

// Methods
bool Teacher::LikesBlue() { return (male || age >= 100); }
