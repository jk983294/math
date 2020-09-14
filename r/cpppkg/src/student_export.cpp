// Include Rcpp system header file (e.g. <>)
#include <Rcpp.h>

// Include our definition of the student file (e.g. "")
#include "student.h"

// std::string GetName(const Student& ps) {
//    return ps.GetName();
//}

// Expose (some of) the Student class
RCPP_MODULE(RcppStudentEx) {
    Rcpp::class_<Student>("Student")
        .constructor<std::string, int, bool>()
        .method("GetName", &Student::GetName)
        .method("GetAge", &Student::GetAge)
        .method("IsMale", &Student::IsMale)
        .method("GetFavoriteNumbers", &Student::GetFavoriteNumbers)
        .method("GetNumbers", &Student::GetNumbers)
        .method("LikesBlue", &Student::LikesBlue);

    Rcpp::class_<BadStudent>("BadStudent")
        .constructor<std::string, int, bool>()
        .method("GetNumbers", &BadStudent::GetNumbers);

    Rcpp::class_<Teacher>("Teacher")
        .constructor<std::string, int, bool>()
        .method("GetName", &Teacher::GetName)
        .method("GetAge", &Teacher::GetAge)
        .method("IsMale", &Teacher::IsMale)
        .method("LikesBlue", &Teacher::LikesBlue);

    // Rcpp::function( "GetName", &GetName );
}
