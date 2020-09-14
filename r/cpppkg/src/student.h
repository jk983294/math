//  Student.h

#ifndef Student_H
#define Student_H

#include <string>
#include <vector>

class Student {
public:
    // Constructor
    Student(std::string name, int age, bool male);

    // Getters
    std::string GetName();
    int GetAge();
    bool IsMale();
    std::vector<int> GetFavoriteNumbers();
    void SetFavoriteNumbers(const std::vector<int>& data);
    virtual std::vector<int> GetNumbers();

    // Methods
    bool LikesBlue();

private:
    // Member variables
    std::string name;
    int age;
    bool male;
    std::vector<int> favoriteNumbers;
};

class BadStudent : public Student {
public:
    BadStudent(std::string name, int age, bool male);
    virtual std::vector<int> GetNumbers();

    using Student::GetName;
};

class Teacher {
public:
    // Constructor
    Teacher(std::string name, int age, bool male);

    // Getters
    std::string GetName();
    int GetAge();
    bool IsMale();

    // Methods
    bool LikesBlue();

private:
    std::string name;
    int age;
    bool male;
};

#endif /* Student_H */
