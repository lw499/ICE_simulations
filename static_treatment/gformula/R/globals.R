if(getRversion() >= "2.15.1"){
  # To remove 'no visible binding for global variable ...' in dplyr and data.table commands
  utils::globalVariables(c('.', 't0', 'id', 'J'))
}
