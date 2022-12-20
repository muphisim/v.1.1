//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file configuration.cpp
  \brief This file contains a function to create the timestamp for the macros
*/
#include <time.h>
#include <configuration.h>

extern char *Timestamp() {
    time_t timestamp = time(NULL);
    char *date = ctime(&timestamp);
    date[strlen(date) - 1] = '\0';
    return date;
};

