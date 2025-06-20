#pragma once

#include <iostream>    // Required for input/output operations (e.g., std::cout)
#include <string>      // Required for std::string
#include <sstream>     // Required for std::stringstream (for efficient string building)
#include <utility>     // Required for std::forward (for perfect forwarding in templates)


namespace rad{
  namespace utils{
    /**
     * @brief Helper struct to append arguments with a comma separator.
     * This is used internally by the createFunctionCallString template.
     * @tparam T The type of the current argument.
     */
    template <typename T>
      struct ArgAppender {
	std::stringstream& ss; // Reference to the stringstream where the string is being built
	bool& firstArg;        // Reference to a boolean flag to track if it's the first argument

	/**
	 * @brief Constructor for ArgAppender.
	 * @param stream The stringstream to append to.
	 * @param isFirst Reference to the first argument flag.
	 */
      ArgAppender(std::stringstream& stream, bool& isFirst) : ss(stream), firstArg(isFirst) {}

	/**
	 * @brief Overloaded operator() to append an argument.
	 * It adds a comma before the argument if it's not the first one.
	 * @param arg The argument to append.
	 */
	void operator()(T&& arg) const {
	  // If it's not the first argument, add a comma and a space for readability.
	  if (!firstArg) {
            ss << ", ";
	  }
	  // Append the argument to the stringstream.
	  ss << std::forward<T>(arg);
	  // After appending, it's no longer the first argument.
	  firstArg = false;
	}
      };

    /**
     * @brief Creates a string representation of a function call.
     * This version handles functions with zero arguments.
     * @param funcName The name of the function.
     * @return A string representing the function call (e.g., "myFunc()").
     */
    std::string createFunctionCallString(const std::string& funcName) {
      return funcName + "()";
    }

    /**
     * @brief Creates a string representation of a function call.
     * This template version handles functions with one or more arguments.
     * @tparam Args The types of the arguments (deduced automatically).
     * @param funcName The name of the function.
     * @param args The arguments to be included in the function call string.
     * These can be any types that can be streamed to std::stringstream
     * (e.g., std::string, int, double, etc.).
     * @return A string representing the function call (e.g., "myFunc(arg1, arg2)").
     */
    template <typename... Args>
      std::string createFunctionCallString(const std::string& funcName, Args&&... args) {
      std::stringstream ss; // Create a stringstream to build the result string.
      ss << funcName << "("; // Start with the function name and opening parenthesis.

      bool firstArg = true; // Flag to manage comma separation for arguments.

      // Use a fold expression (C++17 and later) to iterate over arguments.
      // For each argument, an ArgAppender object is created and invoked,
      // which handles adding commas and appending the argument to the stringstream.
      (ArgAppender<Args>(ss, firstArg)(std::forward<Args>(args)), ...);

      ss << ")"; // End with the closing parenthesis.
      return ss.str(); // Return the built string.
    }

   /**
     * @brief Creates a string representation of a function call.
     * This template version handles a vector of string arguments.
     * @param funcName The name of the function.
     * @param args The arguments to be included in the function call string.
     * These can be any types that can be streamed to std::stringstream
     * (e.g., std::string, int, double, etc.).
     * @return A string representing the function call (e.g., "myFunc(arg1, arg2)").
     */
    std::string createFunctionCallStringFromVec(const std::string& funcName, const  std::vector<std::string>&  args) {
      std::stringstream ss; // Create a stringstream to build the result string.
      ss << funcName << "("<<args[0]; // Start with the function name and opening parenthesis.
      if(args.size()>1){
	std::for_each(args.begin()+1, args.end(),
		       [&ss](const std::string& s) {
			 ss << "," << s;
		       });
      }
      ss << ")"; // End with the closing parenthesis.
      return ss.str(); // Return the built string.
    }

    /**
     * @brief Replaces all occurrences of a specified substring within a string.
     *
     * This function iterates through the input string, finding all instances of
     * 'oldSubstr' and replacing them with 'newSubstr'.
     *
     * @param str The original string in which replacements will be made.
     * @param oldSubstr The substring to be replaced.
     * @param newSubstr The substring to replace 'oldSubstr' with.
     * @return A new string with all occurrences replaced.
     */
    std::string replaceAll(std::string& str, const std::string& oldSubstr, const std::string& newSubstr) {
      // Start searching from the beginning of the string.
      size_t pos = 0;

      // Loop until no more occurrences of oldSubstr are found.
      while ((pos = str.find(oldSubstr, pos)) != std::string::npos) {
        // Replace the found occurrence.
        // str.replace(position, length_of_old_substring, new_substring_content)
        str.replace(pos, oldSubstr.length(), newSubstr);

        // Advance the search position by the length of the new substring.
        // This is crucial to avoid infinite loops if newSubstr contains oldSubstr
        // (e.g., replacing "a" with "aa") and to continue searching
        // after the newly inserted text.
        pos += newSubstr.length();
      }
      // Return the modified string.
      return str;
    }

/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{item1, item2, item3}").
 */
std::string combineVectorToString(const std::vector<std::string>& stringVector) {
    std::stringstream ss; // Create a stringstream for efficient string building.
    ss << "{";            // Prepend with an opening curly brace.

    // Use a loop to iterate through the vector elements.
    for (size_t i = 0; i < stringVector.size(); ++i) {
        ss << stringVector[i]; // Append the current string.

        // If it's not the last element, append a comma and a space.
        if (i < stringVector.size() - 1) {
            ss << ", ";
        }
    }

    ss << "}"; // Append with a closing curly brace.

    return ss.str(); // Return the final combined string.
}


  }
}
