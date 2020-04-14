program: functions? data? tdata? params? tparams? model? generated?;

functions: 'functions' function_decls;
data: 'data' var_decls;
tdata: 'transformed data' var_decls_statements;
params: 'parameters' var_decls;
tparams: 'transformed parameters' var_decls_statements;
model: 'model' var_decls_statements;
generated: 'generated quantities' var_decls_statements;
function_decls: '{' function_decl* '}';
var_decls: '{' var_decl* '}';
var_decls_statements: '{' var_decl* statement* '}';


function_decl: return_type identifier '(' (parameter_decl (',' parameter_decl)*)? ')'
                  statement;

return_type: 'void' | unsized_type;
parameter_decl: 'data'? unsized_type identifier;
unsized_type: basic_type unsized_dims?;
basic_type: 'int' | 'real' | 'vector' | 'row_vector' | 'matrix';
unsized_dims: '['  ','*  ']';




var_decl: var_type variable dims? ('=' expression)? ';';

var_type: 'int' range_constraint
           | 'real' constraint
           | 'vector' constraint '[' expression ']'
           | 'ordered' '[' expression ']'
           | 'positive_ordered' '[' expression ']'
           | 'simplex' '[' expression ']'
           | 'unit_vector' '[' expression ']'
           | 'row_vector' constraint '[' expression ']'
           | 'matrix' constraint '[' expression ',' expression ']'
           | 'cholesky_factor_corr' '[' expression ']'
           | 'cholesky_factor_cov' '[' expression (',' expression)? ']'
           | 'corr_matrix' '[' expression ']'
           | 'cov_matrix' '[' expression ']';

constraint: range_constraint | ('<' offset_multiplier '>');

range_constraint: ('<' range '>')?;

range: 'lower' '=' constr_expression ',' 'upper' '=' constr_expression
        | 'lower' '=' constr_expression
        | 'upper' '=' constr_expression;


offset_multiplier: 'offset' '=' constr_expression ','
                      'multiplier' '=' constr_expression
            | 'offset' '=' constr_expression
            | 'multiplier' '=' constr_expression;

dims: '['  expressions ']';

variable: identifier;

identifier: "[a-zA-Z][a-zA-Z0-9_]*";



expressions: (expression (',' expression)*)?;

expression: expression '?' expression ':' expression
             | expression infixOp expression
             | prefixOp expression
             | expression postfixOp
             | common_expression;

constr_expression: constr_expression arithmeticInfixOp constr_expression
                    | prefixOp constr_expression
                    | constr_expression postfixOp
                    | constr_expression '[' indexes ']'
                    | common_expression;

common_expression : real_literal
    | variable
    | '{' expressions '}'
    | '[' expressions ']'
    | variable '[' expressions ']'
    | function_literal '(' expressions? ')'
    | function_literal '(' expression ('|' (expression (',' expression)*)?)? ')'
    //| 'integrate_1d' '(' function_literal (',' expression)@5:6 ')'
    //| 'integrate_ode' '(' function_literal (',' expression)@6 ')'
    //| 'integrate_ode_rk45' '(' function_literal (',' expression)@6:9 ')'
    //| 'integrate_ode_bdf' '(' function_literal (',' expression)@6:9 ')'
    //| 'algebra_solver' '(' function_literal (',' expression)@4:7 ')'
    //| 'map_rect' '(' function_literal (',' expression)@4 ')'
    | '(' expression ')';

prefixOp: ('!' | '-' | '+' | '^');

postfixOp: '\'';

infixOp: arithmeticInfixOp | logicalInfixOp;

arithmeticInfixOp: ('+' | '-' | '*' | '/' | '%' | '\\' | '.*' | './');

logicalInfixOp: ('||' | '&&' | '==' | '!=' | '<' | '<=' | '>' | '>=');

index: (expression | expression ':' | ':' expression
        | expression ':' expression)?;

indexes: (index (',' index)*)?;

integer_literal: "[0-9]+";

real_literal: integer_literal '.' "[0-9]*" exp_literal?
               | '.' "[0-9]+" exp_literal?
               | integer_literal exp_literal?;

exp_literal: ('e' | 'E') ('+' | '-')? integer_literal;

function_literal: identifier;



statement: atomic_statement | nested_statement;

atomic_statement:  lhs assignment_op expression ';'
   | sampling_statement
   | function_literal '(' expressions ')' ';'
   | 'increment_log_prob' '(' expression ')' ';'
   | 'target' '+=' expression ';'
   | 'break' ';'
   | 'continue' ';'
   | 'print' '(' ((expression | string_literal) (',' (expression | string_literal))*)? ')' ';'
   | 'reject' '(' ((expression | string_literal) (',' (expression | string_literal))*)? ')' ';'
   | 'return' expression ';'
   | ';';

sampling_statement: expression '~' identifier '(' expressions ')' truncation? ';';

assignment_op: '<-' | '=' | '+=' | '-=' | '*=' | '/=' | '.*=' | './=';

//string_literal: '"' char* '"';
string_literal: "\"([^\"\\]|\\[^])*\"";

truncation: 'T' '[' ?expression ',' ?expression ']';

lhs: identifier ('[' indexes ']')*;

nested_statement: 'if' '(' expression ')' statement
    ('else' 'if' '(' expression ')' statement)*
    ('else' statement)?
  | 'while' '(' expression ')' statement
  | 'for' '(' identifier 'in' expression ':' expression ')' statement
  | 'for' '(' identifier 'in' expression ')' statement
  | '{' var_decl* statement+ '}';
